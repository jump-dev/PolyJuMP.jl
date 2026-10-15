export Poly

for op in (:+, :-)
    @eval begin
        function Base.$op(
            a::JuMP.AbstractJuMPScalar,
            p::MP.AbstractPolynomialLike,
        )
            p = p isa MP.AbstractMonomialLike ? p : MB.implicit(p)
            return $op(MP.constant_term(a, p), p)
        end

        function Base.$op(
            p::MP.AbstractPolynomialLike,
            a::JuMP.AbstractJuMPScalar,
        )
            p = p isa MP.AbstractMonomialLike ? p : MB.implicit(p)
            return $op(p, MP.constant_term(a, p))
        end
    end
end

# Resolve the overlap with MP's methods for scalars of the coefficient type.
for Scalar in (JuMP.AbstractJuMPScalar, MOI.AbstractScalarFunction),
    op in (:+, :-)

    @eval begin
        function Base.$op(
            a::T,
            p::Union{MP.AbstractPolynomial{T},MP.AbstractTerm{T}},
        ) where {T<:$Scalar}
            p = MB.implicit(p)
            return $op(MP.constant_term(a, p), p)
        end

        function Base.$op(
            p::Union{MP.AbstractPolynomial{T},MP.AbstractTerm{T}},
            a::T,
        ) where {T<:$Scalar}
            p = MB.implicit(p)
            return $op(p, MP.constant_term(a, p))
        end
    end
end

function Base.:*(a::JuMP.AbstractJuMPScalar, p::MP.AbstractPolynomialLike)
    return MP.left_constant_mult(a, p)
end

function Base.:*(p::MP.AbstractPolynomialLike, a::JuMP.AbstractJuMPScalar)
    return MP.right_constant_mult(p, a)
end

function JuMP.value(f::Function, t::MP.AbstractTerm{<:JuMP.AbstractJuMPScalar})
    return JuMP.value(f, MP.coefficient(t)) * MP.monomial(t)
end

function JuMP.value(t::MP.AbstractTerm{<:JuMP.AbstractJuMPScalar})
    return JuMP.value(JuMP.value, t)
end

function JuMP.value(
    f::Function,
    p::SA.AlgebraElement{<:JuMP.AbstractJuMPScalar},
)
    return SA.AlgebraElement(JuMP.value.(f, SA.coeffs(p)), parent(p))
end

function JuMP.value(p::SA.AlgebraElement{<:JuMP.AbstractJuMPScalar})
    return JuMP.value(JuMP.value, p)
end

abstract type AbstractPoly end

"""
    struct Poly{PB<:AbstractPolynomialBasis} <: AbstractPoly
        polynomial_basis::PB
    end

Polynomial variable ``v^\\top p`` where ``v`` is a vector of new decision variables and ``p`` is a vector of polynomials for the basis `polynomial_basis`.
"""
struct Poly{B<:SA.ExplicitBasis} <: AbstractPoly
    polynomial_basis::B
end
function Poly(x::AbstractVector{<:MP.AbstractPolynomialLike})
    return Poly(MB.SubBasis{MB.Monomial}(x))
end

# Macro
function cvarchecks(_error::Function, info; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    if !isnan(info.start)
        _error(
            "Polynomial variable declaration does not support the start keyword argument.",
        )
    end
    if info.has_fix
        _error(
            "Polynomial variable declaration does not support the form ... == value.",
        )
    end
    if info.has_lb && info.has_ub
        _error(
            "Polynomial variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.",
        )
    end
    if info.has_lb
        _error(
            "Polynomial variable declaration does not support the form ... >= lb.",
        )
    end
    if info.has_ub
        _error(
            "Polynomial variable declaration does not support the form ... <= ub.",
        )
    end
end
function _warnbounds(_error, p::AbstractPoly, info) end
function _warnbounds(_error, p::Poly, info)
    if info.has_lb
        _error(
            "Free polynomial variable declaration does not support the form ... >= 0, use SOSPoly(x) instead of Poly(x) to create Sum of Squares polynomials. Note that SOSPoly(x) creates the polynomial x^T Q x with Q symmetric positive semidefinite while Poly(x) creates the polynomial a^T x so the meaning of the vector of monomial x changes from Poly to SOSPoly.",
        )
    end
end

struct Variable{PT<:AbstractPoly} <: JuMP.AbstractVariable
    p::PT
    binary::Bool
    integer::Bool
end
function JuMP.build_variable(
    _error::Function,
    info::JuMP.VariableInfo,
    p::AbstractPoly;
    extra_kwargs...,
)
    cvarchecks(_error, info; extra_kwargs...)
    _warnbounds(_error, p, info)
    return Variable(p, info.binary, info.integer)
end

# Free polynomial
function JuMP.add_variable(
    model::JuMP.AbstractModel,
    v::Variable{<:Poly},
    name::String = "",
)
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
    return MB.algebra_element(_newvar, v.p.polynomial_basis)
end
