export Poly

function JuMP.resultvalue(t::AbstractTerm{<:JuMP.AbstractJuMPScalar})
    JuMP.resultvalue(coefficient(t)) * monomial(t)
end
function JuMP.resultvalue(p::AbstractPolynomialLike{<:JuMP.AbstractJuMPScalar})
    polynomial(JuMP.resultvalue.(terms(p)), MultivariatePolynomials.SortedUniqState())
end

abstract type AbstractPoly end

"""
    struct Poly{PB<:AbstractPolynomialBasis} <: AbstractPoly
        polynomial_basis::PB
    end

Polynomial variable ``v^\\top p`` where ``v`` is a vector of new decision variables and ``p`` is a vector of polynomials for the basis `polynomial_basis`.
"""
struct Poly{PB<:AbstractPolynomialBasis} <: AbstractPoly
    polynomial_basis::PB
end
Poly(x::AbstractVector{<:MultivariatePolynomials.AbstractPolynomialLike}) = Poly(MonomialBasis(x))

# Macro
function cvarchecks(_error::Function, info; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    if !isnan(info.start)
        _error("Polynomial variable declaration does not support the start keyword argument.")
    end
    if info.hasfix
        _error("Polynomial variable declaration does not support the form ... == value.")
    end
    if info.haslb && info.hasub
        _error("Polynomial variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.")
    end
    if info.haslb
        _error("Polynomial variable declaration does not support the form ... >= lb.")
    end
    if info.hasub
        _error("Polynomial variable declaration does not support the form ... <= ub.")
    end
end
function _warnbounds(_error, p::AbstractPoly, info) end
function _warnbounds(_error, p::Poly, info)
    if info.haslb
        _error("Free polynomial variable declaration does not support the form ... >= 0, use SOSPoly(x) instead of Poly(x) to create Sum of Squares polynomials. Note that SOSPoly(x) creates the polynomial x^T Q x with Q symmetric positive semidefinite while Poly(x) creates the polynomial a^T x so the meaning of the vector of monomial x changes from Poly to SOSPoly.")
    end
end

struct Variable{PT<:AbstractPoly} <: JuMP.AbstractVariable
    p::PT
    binary::Bool
    integer::Bool
end
function JuMP.buildvariable(_error::Function, info::JuMP.VariableInfo, p::AbstractPoly; extra_kwargs...)
    cvarchecks(_error, info; extra_kwargs...)
    _warnbounds(_error, p, info)
    Variable(p, info.binary, info.integer)
end
function JuMP.addvariable(m::JuMP.AbstractModel, v::Variable, name::String)
    createpoly(m, getdefault(m, v.p), v.binary, v.integer)
end
