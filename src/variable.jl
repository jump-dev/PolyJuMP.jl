export Poly

function JuMP.value(
    t::AbstractTerm{<:JuMP.AbstractJuMPScalar},
    f::Function = JuMP.value,
)
    return JuMP.value(MultivariatePolynomials.coefficient(t), f) * monomial(t)
end
function JuMP.value(
    p::AbstractPolynomialLike{<:JuMP.AbstractJuMPScalar},
    f::Function = JuMP.value,
)
    return polynomial(
        JuMP.value.(terms(p), f),
        MultivariatePolynomials.SortedUniqState(),
    )
end

abstract type AbstractPoly end

"""
    struct Poly{PB<:AbstractPolynomialBasis} <: AbstractPoly
        polynomial_basis::PB
    end

Polynomial variable ``v^\\top p`` where ``v`` is a vector of new decision variables and ``p`` is a vector of polynomials for the basis `polynomial_basis`.
"""
struct Poly{PB<:MB.AbstractPolynomialBasis} <: AbstractPoly
    polynomial_basis::PB
end
Poly(x::AbstractVector{<:MultivariatePolynomials.AbstractPolynomialLike}) = Poly(MB.MonomialBasis(x))

# Macro
function cvarchecks(_error::Function, info; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    if !isnan(info.start)
        _error("Polynomial variable declaration does not support the start keyword argument.")
    end
    if info.has_fix
        _error("Polynomial variable declaration does not support the form ... == value.")
    end
    if info.has_lb && info.has_ub
        _error("Polynomial variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.")
    end
    if info.has_lb
        _error("Polynomial variable declaration does not support the form ... >= lb.")
    end
    if info.has_ub
        _error("Polynomial variable declaration does not support the form ... <= ub.")
    end
end
function _warnbounds(_error, p::AbstractPoly, info) end
function _warnbounds(_error, p::Poly, info)
    if info.has_lb
        _error("Free polynomial variable declaration does not support the form ... >= 0, use SOSPoly(x) instead of Poly(x) to create Sum of Squares polynomials. Note that SOSPoly(x) creates the polynomial x^T Q x with Q symmetric positive semidefinite while Poly(x) creates the polynomial a^T x so the meaning of the vector of monomial x changes from Poly to SOSPoly.")
    end
end

struct Variable{PT<:AbstractPoly} <: JuMP.AbstractVariable
    p::PT
    binary::Bool
    integer::Bool
end
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, p::AbstractPoly; extra_kwargs...)
    cvarchecks(_error, info; extra_kwargs...)
    _warnbounds(_error, p, info)
    return Variable(p, info.binary, info.integer)
end

# Free polynomial
function JuMP.add_variable(model::JuMP.AbstractModel, v::Variable{<:Poly},
                           name::String="")
    function _newvar(i)
        vref = VariableRef(model)
        if v.binary
            JuMP.set_binary(vref)
        end
        if v.integer
            JuMP.set_integer(vref)
        end
        return vref
    end
    return polynomial(_newvar, v.p.polynomial_basis)
end
