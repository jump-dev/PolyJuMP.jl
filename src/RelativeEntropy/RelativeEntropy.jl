module RelativeEntropy

import MutableArithmetics as MA
import MultivariatePolynomials as MP
import MathOptInterface as MOI
import JuMP
import PolyJuMP

abstract type AbstractAGECone <: MOI.AbstractVectorSet end

"""
    struct SignomialSAGECone <: MOI.AbstractVectorSet
        α::Matrix{Int}
    end

**S**ums of **A**M/**G**M **E**xponential for signomials.
"""
struct SignomialSAGECone <: AbstractAGECone
    α::Matrix{Int}
end

"""
    struct PolynomialSAGECone <: MOI.AbstractVectorSet
        α::Matrix{Int}
    end

**S**ums of **A**M/**G**M **E**xponential for polynomials.
"""
struct PolynomialSAGECone <: AbstractAGECone
    α::Matrix{Int}
end

struct SignomialAGECone <: AbstractAGECone
    α::Matrix{Int}
    k::Int
end

struct PolynomialAGECone <: AbstractAGECone
    α::Matrix{Int}
    k::Int
end

MOI.dimension(set::AbstractAGECone) = size(set.α, 1)
Base.copy(set::AbstractAGECone) = set

function _exponents_matrix(monos)
    α = Matrix{Int}(undef, length(monos), MP.nvariables(monos))
    for (i, mono) in enumerate(monos)
        exp = MP.exponents(mono)
        for j in eachindex(exp)
            α[i, j] = exp[j]
        end
    end
    return α
end

struct SignomialSAGESet <: PolyJuMP.PolynomialSet end
function JuMP.reshape_set(::SignomialSAGECone, ::PolyJuMP.PolynomialShape)
    return SignomialSAGESet()
end
function JuMP.moi_set(::SignomialSAGESet, monos)
    return SignomialSAGECone(_exponents_matrix(monos))
end

struct PolynomialSAGESet <: PolyJuMP.PolynomialSet end
function JuMP.reshape_set(::PolynomialSAGECone, ::PolyJuMP.PolynomialShape)
    return PolynomialSAGESet()
end
function JuMP.moi_set(::PolynomialSAGESet, monos)
    return PolynomialSAGECone(_exponents_matrix(monos))
end

struct SignomialAGESet{MT<:MP.AbstractMonomial} <: PolyJuMP.PolynomialSet
    monomial::MT
end
function JuMP.reshape_set(
    set::SignomialAGECone,
    shape::PolyJuMP.PolynomialShape,
)
    return SignomialAGESet(shape.monomials[set.k])
end
function JuMP.moi_set(set::SignomialAGESet, monos)
    k = findfirst(isequal(set.monomial), monos)
    return SignomialAGECone(_exponents_matrix(monos), k)
end

struct PolynomialAGESet{MT<:MP.AbstractMonomial} <: PolyJuMP.PolynomialSet
    monomial::MT
end
function JuMP.reshape_set(
    set::PolynomialAGECone,
    shape::PolyJuMP.PolynomialShape,
)
    return PolynomialAGESet(shape.monomials[set.k])
end
function JuMP.moi_set(set::PolynomialAGESet, monos)
    k = findfirst(isequal(set.monomial), monos)
    return PolynomialAGECone(_exponents_matrix(monos), k)
end

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, PolynomialSAGESet)
    return
end

function JuMP.build_constraint(
    _error::Function,
    p,
    set::Union{
        SignomialSAGESet,
        PolynomialSAGESet,
        SignomialAGESet,
        PolynomialAGESet,
    };
    kws...,
)
    coefs = PolyJuMP.non_constant_coefficients(p)
    monos = MP.monomials(p)
    cone = JuMP.moi_set(set, monos)
    shape = PolyJuMP.PolynomialShape(monos)
    return PolyJuMP.bridgeable(
        JuMP.VectorConstraint(coefs, cone, shape),
        JuMP.moi_function_type(typeof(coefs)),
        typeof(cone),
    )
end

const APL{T} = MP.AbstractPolynomialLike{T}

"""
    struct Decomposition{T, PT}

Represents a SAGE decomposition.
"""
struct Decomposition{T,P<:APL{T}} <: APL{T}
    polynomials::Vector{P}
    function Decomposition(ps::Vector{P}) where {T,P<:APL{T}}
        return new{T,P}(ps)
    end
end

function Base.show(io::IO, p::Decomposition)
    for (i, q) in enumerate(p.polynomials)
        print(io, "(")
        print(io, q)
        print(io, ")")
        if i != length(p.polynomials)
            print(io, " + ")
        end
    end
end

function MP.polynomial(d::Decomposition)
    return sum(d.polynomials)
end

"""
    struct DecompositionAttribute{T} <: MOI.AbstractConstraintAttribute
        tol::T
    end

A constraint attribute for the [`Decomposition`](@ref) of a constraint.
"""
struct DecompositionAttribute{T} <: MOI.AbstractConstraintAttribute
    tol::T
    result_index::Int
end
function DecompositionAttribute(tol::Real)
    return DecompositionAttribute(tol, 1)
end

function decomposition(
    con_ref::JuMP.ConstraintRef;
    tol::Real,
    result_index::Int = 1,
)
    monos = con_ref.shape.monomials
    attr = DecompositionAttribute(tol, result_index)
    return Decomposition([
        MP.polynomial(a, monos) for
        a in MOI.get(JuMP.owner_model(con_ref), attr, con_ref)
    ])
end

MOI.is_set_by_optimize(::DecompositionAttribute) = true

include("bridges/sage.jl")

function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SignomialSAGECone},
)
    return [(SAGEBridge, PolyJuMP._coef_type(F))]
end

include("bridges/age.jl")

function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SignomialAGECone},
)
    return [(AGEBridge, PolyJuMP._coef_type(F))]
end

include("bridges/signomial.jl")

function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:Union{PolynomialSAGECone,PolynomialAGECone}},
)
    return [(SignomialBridge, PolyJuMP._coef_type(F))]
end

end
