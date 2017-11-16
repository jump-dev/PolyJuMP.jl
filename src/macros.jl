using JuMP
import JuMP: validmodel, addtoexpr_reorder
using Base.Meta

export Poly, @set

function JuMP.resultvalue(t::AbstractTerm{<:JuMP.AbstractJuMPScalar})
    JuMP.resultvalue(coefficient(t)) * monomial(t)
end
function JuMP.resultvalue(p::AbstractPolynomialLike{<:JuMP.AbstractJuMPScalar})
    polynomial(JuMP.resultvalue.(terms(p)), MultivariatePolynomials.SortedUniqState())
end

abstract type AbstractPoly end

# x is a vector of monomials to be used to construct a polynomial variable
# if MT is Gram, x represents the monomials of the form x^T Q x
# if MT is Classic, it represents the monomials of the form a^T x
# if MT is Default, it depends on whether the polynomials is constructed as nonnegative or not:
# For a nonnegative polynomial, it corresponds to Gram, otherwise it corresponds to Classic.
struct Poly{P, MT, MV} <: AbstractPoly
    x::MV
end
Poly{P, MT}(x::MV) where {P, MT, MV} = Poly{P, MT, MV}(x)
Poly{P}(x) where P = Poly{P, :Default}(x)
Poly(x) = Poly{false}(x)

function JuMP.variabletype(m::Model, p::AbstractPoly)
    getpolymodule(m).polytype(m, p)
end
function cvarchecks(_error::Function, haslb::Bool, lowerbound::Number, hasub::Bool, hasfix::Bool, hasstart::Bool; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    if hasfix || hasstart
        _error("Polynomial variable declaration does not support the form ... == value.")
    end
    if haslb && hasub
        _error("Polynomial variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.")
    end
    if haslb && !iszero(lowerbound)
        _error("Polynomial variable declaration does not support the form ... >= lb with nonzero lb.")
    end
    if hasub
        _error("Polynomial variable declaration does not support the form ... <= ub.")
    end
end
function _warnbounds(_error, p::AbstractPoly, haslb, hasub) end
function _warnbounds(_error, p::Poly{false}, haslb, hasub)
    if haslb
        _error("Free polynomial variable declaration does not support the form ... >= 0, use SOSPoly(x) instead of Poly(x) to create Sum of Squares polynomials. Note that SOSPoly(x) creates the polynomial x^T Q x with Q symmetric positive semidefinite while Poly(x) creates the polynomial a^T x so the meaning of the vector of monomial x changes from Poly to SOSPoly.")
    end
end
function JuMP.constructvariable!(m::Model, p::AbstractPoly, _error::Function,
                                 haslb::Bool, lowerbound::Number,
                                 hasub::Bool, upperbound::Number,
                                 hasfix::Bool, fixedvalue::Number,
                                 binary::Bool, integer::Bool, name::AbstractString,
                                 hasstart::Bool, start::Number; extra_kwargs...)
    cvarchecks(_error, haslb, lowerbound, hasub, hasfix, hasstart; extra_kwargs...)
    _warnbounds(_error, p, haslb, hasub)
    getpolymodule(m).createpoly(m, p, binary, integer)
end

using MathOptInterface
const MOI = MathOptInterface

function JuMP.constructconstraint!(p, s::MOI.EqualTo)
    PolyConstraint(p-s.value, ZeroPoly())
end
function JuMP.constructconstraint!(p, s::MOI.GreaterThan)
    PolyConstraint(p-s.lower, NonNegPoly())
end
function JuMP.constructconstraint!(p, s::MOI.LessThan)
    PolyConstraint(s.upper-p, NonNegPoly())
end

function JuMP.constructconstraint!{PolyT<:AbstractPolynomialLike}(p::Union{PolyT, AbstractMatrix{PolyT}}, s)
    PolyConstraint(p, s)
end
# there is already a method for AbstractMatrix in PSDCone in JuMP so we need a more specific here to avoid ambiguity
function JuMP.constructconstraint!{PolyT<:AbstractPolynomialLike}(p::AbstractMatrix{PolyT}, s::PSDCone)
    PolyConstraint(p, s)
end
