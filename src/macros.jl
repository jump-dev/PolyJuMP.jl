export Poly

function JuMP.getvalue(t::AbstractTerm{<:JuMP.AbstractJuMPScalar})
    getvalue(coefficient(t)) * monomial(t)
end
function JuMP.getvalue(p::AbstractPolynomialLike{<:JuMP.AbstractJuMPScalar})
    polynomial(getvalue.(terms(p)), MultivariatePolynomials.SortedUniqState())
end

abstract type AbstractPoly end

# x is a vector of monomials to be used to construct a polynomial variable
# if MS is Gram, x represents the monomials of the form x^T Q x
# if MS is Classic, it represents the monomials of the form a^T x
# if MS is Default, it depends on whether the polynomials is constructed as nonnegative or not:
# For a nonnegative polynomial, it corresponds to Gram, otherwise it corresponds to Classic.
struct Poly{P, MS, MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPoly
    x::MV
end
Poly{P, MS}(x::AbstractVector{MT}) where {P, MS, MT<:MultivariatePolynomials.AbstractMonomial} = Poly{P, MS, MT, typeof(x)}(x)
Poly{P, MS}(x) where {P, MS} = Poly{P, MS}(monovec(x))
Poly{P}(x) where P = Poly{P, :Default}(x)
Poly(x) = Poly{false}(x)

JuMP.variabletype(m::JuMP.Model, p::Poly{true}) = JuMP.variabletype(m, getdefault(m, p))

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
    if lowerbound != -Inf && lowerbound != 0
        _error("Polynomial variable declaration does not support the form ... >= lb with nonzero lb.")
    end
    if upperbound != Inf
        _error("Polynomial variable declaration does not support the form ... <= ub.")
    end
end
function _warnbounds(_error, p::AbstractPoly, lowerbound, upperbound) end
function _warnbounds(_error, p::Poly{false}, lowerbound, upperbound)
    if lowerbound != -Inf
        _error("Free polynomial variable declaration does not support the form ... >= 0, use SOSPoly(x) instead of Poly(x) to create Sum of Squares polynomials. Note that SOSPoly(x) creates the polynomial x^T Q x with Q symmetric positive semidefinite while Poly(x) creates the polynomial a^T x so the meaning of the vector of monomial x changes from Poly to SOSPoly.")
    end
end
function JuMP.constructvariable!(m::Model, p::AbstractPoly, _error::Function, lowerbound::Number, upperbound::Number, category::Symbol, basename::AbstractString, start::Number; extra_kwargs...)
    cvarchecks(_error, lowerbound, upperbound, start; extra_kwargs...)
    _warnbounds(_error, p, lowerbound, upperbound)
    createpoly(m, getdefault(m, p), category == :Default ? :Cont : category)
end

function JuMP.constructconstraint!(p::AbstractPolynomialLike, sense::Symbol)
    JuMP.constructconstraint!(sense == :(<=) ? -p : p, sense == :(==) ? ZeroPoly() : NonNegPoly())
end

function JuMP.constructconstraint!(p::Union{AbstractPolynomialLike, AbstractMatrix{<:AbstractPolynomialLike}}, s)
    PolyConstraint(p, s)
end
# there is already a method for AbstractMatrix in PSDCone in JuMP so we need a more specific here to avoid ambiguity
function JuMP.constructconstraint!(p::AbstractMatrix{<:AbstractPolynomialLike}, s::PSDCone)
    PolyConstraint(p, NonNegPolyMatrix())
end
