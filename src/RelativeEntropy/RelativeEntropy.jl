module RelativeEntropy

import MutableArithmetics as MA
import MultivariatePolynomials as MP
import MathOptInterface as MOI
import JuMP
import PolyJuMP

"""
    struct SAGECone <: MOI.AbstractVectorSet
        α::Matrix{Int}
    end

**S**ums of **A**M/**G**M **E**xponential.
"""
struct SAGECone <: MOI.AbstractVectorSet
    α::Matrix{Int}
end
MOI.dimension(set::SAGECone) = size(set.α, 1)
Base.copy(set::SAGECone) = set

struct SAGESet <: PolyJuMP.PolynomialSet end
JuMP.reshape_set(::SAGECone, ::PolyJuMP.PolynomialShape) = SAGESet()
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
JuMP.moi_set(::SAGESet, monos) = SAGECone(_exponents_matrix(monos))

struct AGECone <: MOI.AbstractVectorSet
    α::Matrix{Int}
    k::Int
end
MOI.dimension(set::AGECone) = size(set.α, 1)
Base.copy(set::AGECone) = set

struct AGESet{MT<:MP.AbstractMonomial} <: PolyJuMP.PolynomialSet
    monomial::MT
end
function JuMP.reshape_set(set::AGECone, shape::PolyJuMP.PolynomialShape)
    return AGESet(shape.monomials[set.k])
end
function JuMP.moi_set(set::AGESet, monos)
    k = findfirst(isequal(set.monomial), monos)
    return AGECone(_exponents_matrix(monos), k)
end

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, SAGESet)
    return
end

function JuMP.build_constraint(
    _error::Function,
    p,
    set::Union{SAGESet,AGESet};
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

include("bridges/sage.jl")

function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SAGECone},
)
    return [(SAGEBridge, PolyJuMP._coef_type(F))]
end

include("bridges/age.jl")

function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:AGECone},
)
    return [(AGEBridge, PolyJuMP._coef_type(F))]
end

end
