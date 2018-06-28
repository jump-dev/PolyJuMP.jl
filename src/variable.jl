export Poly

function JuMP.getvalue(t::AbstractTerm{<:JuMP.AbstractJuMPScalar})
    getvalue(coefficient(t)) * monomial(t)
end
function JuMP.getvalue(p::AbstractPolynomialLike{<:JuMP.AbstractJuMPScalar})
    polynomial(getvalue.(terms(p)), MultivariatePolynomials.SortedUniqState())
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
function cvarchecks(_error::Function, lowerbound::Number, upperbound::Number, start::Number; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    if !isnan(start)
        _error("Polynomial variable declaration does not support the start keyword argument.")
    end
    if lowerbound != -Inf && upperbound != Inf
        if lowerbound == upperbound
            _error("Polynomial variable declaration does not support the form ... == value.")
        else
            _error("Polynomial variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.")
        end
    end
    if lowerbound != -Inf
        _error("Polynomial variable declaration does not support the form ... >= lb.")
    end
    if upperbound != Inf
        _error("Polynomial variable declaration does not support the form ... <= ub.")
    end
end
function _warnbounds(_error, p::AbstractPoly, lowerbound, upperbound) end
function _warnbounds(_error, p::Poly, lowerbound, upperbound)
    if lowerbound != -Inf
        _error("Free polynomial variable declaration does not support the form ... >= 0, use SOSPoly(x) instead of Poly(x) to create Sum of Squares polynomials. Note that SOSPoly(x) creates the polynomial x^T Q x with Q symmetric positive semidefinite while Poly(x) creates the polynomial a^T x so the meaning of the vector of monomial x changes from Poly to SOSPoly.")
    end
end
function JuMP.constructvariable!(m::Model, p::AbstractPoly, _error::Function, lowerbound::Number, upperbound::Number, category::Symbol, basename::AbstractString, start::Number; extra_kwargs...)
    cvarchecks(_error, lowerbound, upperbound, start; extra_kwargs...)
    _warnbounds(_error, p, lowerbound, upperbound)
    createpoly(m, getdefault(m, p), category == :Default ? :Cont : category)
end
