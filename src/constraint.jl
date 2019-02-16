abstract type PolynomialSet end

function JuMP.in_set_string(print_mode, set::PolynomialSet)
    return string(JuMP._math_symbol(print_mode, :in), ' ', set)
end

struct ZeroPoly <: PolynomialSet end
struct NonNegPoly <: PolynomialSet end
struct PosDefPolyMatrix <: PolynomialSet end

function JuMP.function_string(::Type{JuMP.REPLMode},
                              p::MultivariatePolynomials.APL)
    return sprint(show, MIME"text/plain"(), p)
end
function JuMP.function_string(::Type{JuMP.IJuliaMode},
                              p::MultivariatePolynomials.APL)
    # `show` prints `$$` around what `_show` prints.
    return sprint(MultivariatePolynomials._show, MIME"text/latex"(), p)
end

### Shapes for polynomial/moments primal-dual pair ###

# Inspired from `JuMP.dual_shape` docstring example
struct PolynomialShape{MT <: AbstractMonomial,
                       MVT <: AbstractVector{MT}} <: JuMP.AbstractShape
    monomials::MVT
end
function JuMP.reshape_vector(x::Vector, shape::PolynomialShape)
    return polynomial(x, shape.monomials)
end
struct MomentsShape{MT <: AbstractMonomial,
                    MVT <: AbstractVector{MT}} <: JuMP.AbstractShape
    monomials::MVT
end
function JuMP.reshape_vector(x::Vector, shape::MomentsShape)
    return measure(x, shape.monomials)
end
JuMP.dual_shape(shape::PolynomialShape) = MomentsShape(shape.monomials)
JuMP.dual_shape(shape::MomentsShape) = PolynomialShape(shape.monomials)

JuMP.reshape_set(::ZeroPolynomialSet, ::PolynomialShape) = ZeroPoly()
function JuMP.moi_set(::ZeroPoly, monos::AbstractVector{<:AbstractMonomial};
                      domain::AbstractSemialgebraicSet=FullSpace(),
                      basis=MonomialBasis)
    return ZeroPolynomialSet(domain, basis, monos)
end

### @constraint/@SDconstraint macros ###

non_constant(a::Vector{<:Number}) = convert(Vector{AffExpr}, a)
non_constant(a::Vector{<:JuMP.AbstractJuMPScalar}) = a
non_constant_coefficients(p) = non_constant(coefficients(p))

## ZeroPoly
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::ZeroPoly;
                               domain::AbstractSemialgebraicSet=FullSpace(),
                               kws...)
    coefs = non_constant_coefficients(p)
    monos = monomials(p)
    if domain isa BasicSemialgebraicSet
        # p(x) = 0 for all x in a basic semialgebraic set. We replace it by
        # p(x) ≤ 0 and p(x) ≥ 0 for all x in the basic semialgebraic set.
        # We need to determine the cone two use for `NonNegPoly` which is stored
        # in the `ext` field of the `model` so we need the `model` hence we need
        # to go to `JuMP.add_constraint`.

        # We need to recreate the full list of keyword arguments. `kws` is a
        # `Iterators.Pairs, `kws.data` is a `NamedTuple` and `kws.itr` is a
        # `Tuple`.
        all_kws = Iterators.Pairs(merge((domain=domain,), kws.data),
                                  (:domain, kws.itr...))
        return Constraint(_error, p, s, all_kws)
    else
        set = JuMP.moi_set(s, monos; domain=domain, kws...)
        constraint = JuMP.VectorConstraint(coefs, set, PolynomialShape(monos))
        bridgeable = BridgeableConstraint(constraint,
                                          PolyJuMP.ZeroPolynomialBridge)
        if !(domain isa FullSpace)
            bridgeable = BridgeableConstraint(
                constraint, PolyJuMP.ZeroPolynomialInAlgebraicSetBridge)
        end
        return bridgeable
    end
end
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.EqualTo; kws...)
    return JuMP.build_constraint(_error, p-s.value, ZeroPoly(); kws...)
end

## NonNegPoly and PosDefPolyMatrix
# The `model` is not given in `JuMP.build_constraint` so we create a custom
# `Constraint` object and transform the `set` in `JuMP.add_constraint`.
struct Constraint{FT, PT, ST<:PolynomialSet, KWT<:Iterators.Pairs} <: JuMP.AbstractConstraint
    _error::FT
    polynomial_or_matrix::PT
    set::ST
    kws::KWT
end
function JuMP.build_constraint(_error::Function, polynomial_or_matrix,
                               set::Union{NonNegPoly, PosDefPolyMatrix};
                               kws...)
    return Constraint(_error, polynomial_or_matrix, set, kws)
end

# FIXME the domain will not appear in the printing, it should be a field of
#       `ZeroPoly` which is `FullSpace` by default
JuMP.reshape_set(::PlusMinusSet, ::PolynomialShape) = ZeroPoly()
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::Constraint{<:Any, <:Any, ZeroPoly},
                             name::String = "")
    cone = getdefault(model, NonNegPoly())
    coefs = non_constant_coefficients(constraint.polynomial_or_matrix)
    monos = monomials(constraint.polynomial_or_matrix)
    set = PlusMinusSet(JuMP.moi_set(cone, monos; constraint.kws...))
    new_constraint = JuMP.VectorConstraint(coefs, set, PolynomialShape(monos))
    bridgeable = BridgeableConstraint(new_constraint, PolyJuMP.PlusMinusBridge)
    return JuMP.add_constraint(model, new_constraint, name)
end

function JuMP.add_constraint(model::JuMP.Model, constraint::Constraint,
                             name::String = "")
    set = getdefault(model, constraint.set)
    new_constraint = JuMP.build_constraint(constraint._error,
                                           constraint.polynomial_or_matrix, set;
                                           constraint.kws...)
    return JuMP.add_constraint(model, new_constraint, name)
end

# NonNegPoly
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.GreaterThan; kws...)
    return JuMP.build_constraint(_error, p-s.lower, NonNegPoly(); kws...)
end
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.LessThan; kws...)
    return JuMP.build_constraint(_error, s.upper-p, NonNegPoly(); kws...)
end

# PosDefPolyMatrix
# there is already a method for AbstractMatrix in PSDCone in JuMP so we need a
# more specific here to avoid ambiguity
function JuMP.build_constraint(_error::Function,
                               p::AbstractMatrix{<:AbstractPolynomialLike},
                               s::PSDCone; kws...)
    return JuMP.build_constraint(_error, p, PosDefPolyMatrix(); kws...)
end
