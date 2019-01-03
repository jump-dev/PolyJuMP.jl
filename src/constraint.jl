abstract type PolynomialSet end
struct ZeroPoly <: PolynomialSet end
struct NonNegPoly <: PolynomialSet end
struct PosDefPolyMatrix <: PolynomialSet end

#struct Constraint{PT, ST<:PolynomialSet,
#                  BT <: AbstractPolynomialBasis,
#                  DT <: AbstractSemialgebraicSet} <: JuMP.AbstractConstraint
#    p::PT # typically either be a polynomial or a Matrix of polynomials
#    set::ST
#    basis::BT
#    domain::DT
#end
#
## Responsible for getting slack and dual values
#abstract type ConstraintDelegate end
#
#const ConstraintRef{CD<:ConstraintDelegate} = JuMP.ConstraintRef{Model, CD}
#
#Base.show(io::IO, cref::ConstraintRef) = print(io, "PolyJuMP constraint")

### Shapes for polynomial/moments primal-dual pair ###

# Inspired from `JuMP.dual_shape` docstring example
struct PolynomialShape{MT <: AbstractMonomial,
                       MVT <: AbstractVector{MT}} <: JuMP.AbstractShape
    monomials::MVT
end
JuMP.reshape(x::Vector, shape::PolynomialShape) = polynomial(x, shape.monomials)
struct MomentsShape{MT <: AbstractMonomial,
                    MVT <: AbstractVector{MT}} <: JuMP.AbstractShape
    monomials::MVT
end
JuMP.reshape(x::Vector, shape::MomentsShape) = measure(x, shape.monomials)
dual_shape(shape::PolynomialShape) = MomentsShape(shape.monomials)
dual_shape(shape::MomentsShape) = PolynomialShape(shape.monomials)

#getdelegate(c::ConstraintRef) = c.index
#getslack(c::ConstraintRef) = getslack(getdelegate(c))
#JuMP.dual(c::ConstraintRef) = JuMP.dual(getdelegate(c))

### @constraint/@SDconstraint macros ###

## ZeroPoly
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::ZeroPoly;
                               domain::AbstractSemialgebraicSet=FullSpace(),
                               basis=MonomialBasis)
    set = ZeroPolynomialSet(domain, basis, monomials(p))
    return JuMP.build_constraint(_error, coefficients(p), set)
end
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.EqualTo; kws...)
    return JuMP.build_constraint(_error, p-s.value, ZeroPoly(); kws...)
end

## NonNegPoly and PosDefPolyMatrix
# The `model` is not given in `JuMP.build_constraint` so we create a custom
# `Constraint` object and transform the `set` in `JuMP.add_constraint`.
struct Constraint{FT, PT, ST<:PolynomialSet, KWT} <: JuMP.AbstractConstraint
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
