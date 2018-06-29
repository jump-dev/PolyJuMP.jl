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
function cvarchecks(_error::Function, info; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    if info.hasfix || info.hasstart
        _error("Polynomial variable declaration does not support the form ... == value.")
    end
    if info.haslb && info.hasub
        _error("Polynomial variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.")
    end
    if info.haslb && !iszero(info.lowerbound)
        _error("Polynomial variable declaration does not support the form ... >= lb with nonzero lb.")
    end
    if info.hasub
        _error("Polynomial variable declaration does not support the form ... <= ub.")
    end
end
function _warnbounds(_error, p::AbstractPoly, info) end
function _warnbounds(_error, p::Poly{false}, info)
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
    getpolymodule(m).createpoly(m, v.p, v.binary, v.integer)
end

using MathOptInterface
const MOI = MathOptInterface

function JuMP.buildconstraint(_error::Function, p::AbstractPolynomialLike, s::MOI.EqualTo)
    PolyConstraint(p-s.value, ZeroPoly())
end
function JuMP.buildconstraint(_error::Function, p::AbstractPolynomialLike, s::MOI.GreaterThan)
    PolyConstraint(p-s.lower, NonNegPoly())
end
function JuMP.buildconstraint(_error::Function, p::AbstractPolynomialLike, s::MOI.LessThan)
    PolyConstraint(s.upper-p, NonNegPoly())
end

function JuMP.buildconstraint(_error::Function, p::Union{AbstractPolynomialLike, AbstractMatrix{<:AbstractPolynomialLike}}, s)
    PolyConstraint(p, s)
end
# there is already a method for AbstractMatrix in PSDCone in JuMP so we need a more specific here to avoid ambiguity
function JuMP.buildconstraint(_error::Function, p::AbstractMatrix{<:AbstractPolynomialLike}, s::PSDCone)
    PolyConstraint(p, s)
end
