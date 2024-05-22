module SAGE

import MutableArithmetics as MA
import MultivariateBases as MB
import MultivariatePolynomials as MP
import MathOptInterface as MOI
import JuMP
import PolyJuMP

struct Cone{C} <: MOI.AbstractVectorSet
    cone::C
    α::Matrix{Int}
end

function JuMP.reshape_set(c::Cone, ::PolyJuMP.PolynomialShape)
    return c.cone
end

MOI.dimension(set::Cone) = size(set.α, 1)
Base.copy(set::Cone) = set

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

"""
    struct Signomials{M<:Union{Nothing,Int,MP.AbstractMonomial}} <: PolyJuMP.PolynomialSet

**S**ums of **A**M/**G**M **E**xponential for signomials.
"""
struct Signomials{M<:Union{Nothing,Int,MP.AbstractMonomial}} <:
       PolyJuMP.PolynomialSet
    monomial::M
end
Signomials() = Signomials(nothing)
_index(_, ::Nothing) = nothing
_index(basis, mono::MP.AbstractMonomial) = MB.monomial_index(basis, mono)::Int
function JuMP.moi_set(c::Signomials, basis::MB.SubBasis{MB.Monomial})
    monos = basis.monomials
    return Cone(Signomials(_index(basis, c.monomial)), _exponents_matrix(monos))
end

"""
    struct Polynomials{M<:Union{Nothing,Int,MP.AbstractMonomial}} <: PolyJuMP.PolynomialSet

**S**ums of **A**M/**G**M **E**xponential for polynomials.
"""
struct Polynomials{M<:Union{Nothing,Int,MP.AbstractMonomial}} <:
       PolyJuMP.PolynomialSet
    monomial::M
end
Polynomials() = Polynomials(nothing)
function JuMP.moi_set(c::Polynomials, basis::MB.SubBasis{MB.Monomial})
    monos = basis.monomials
    return Cone(
        Polynomials(_index(basis, c.monomial)),
        _exponents_matrix(monos),
    )
end

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, Polynomials)
    return
end

function JuMP.build_constraint(
    _error::Function,
    p,
    set::Union{Signomials,Polynomials};
    kws...,
)
    for (key, _) in kws
        _error("unsupported keyword argument `$key`.")
    end
    coefs = PolyJuMP.non_constant_coefficients(p)
    basis = MB.SubBasis{MB.Monomial}(MP.monomials(p))
    cone = JuMP.moi_set(set, basis)
    shape = PolyJuMP.PolynomialShape(basis)
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
    monos = con_ref.shape.basis.monomials
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
    ::Type{Cone{Signomials{Nothing}}},
)
    return [(SAGEBridge, PolyJuMP._coef_type(F))]
end

include("bridges/age.jl")

function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{Cone{Signomials{Int}}},
)
    return [(AGEBridge, PolyJuMP._coef_type(F))]
end

include("bridges/signomial.jl")

function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{Cone{Polynomials{M}}},
) where {M}
    return [(SignomialsBridge, PolyJuMP._coef_type(F))]
end

end
