export ZeroPoly, NonNegPoly, NonNegPolyMatrix
export getslack

abstract type PolynomialSet end
struct ZeroPoly <: PolynomialSet end
struct NonNegPoly <: PolynomialSet end
struct NonNegPolyMatrix <: PolynomialSet end

struct PolyConstraint{PT, ST<:PolynomialSet} <: JuMP.AbstractConstraint
    p::PT # typically either be a polynomial or a Matrix of polynomials
    set::ST
end

# Responsible for getting slack and dual values
abstract type ConstraintDelegate end
const PolyConstraintRef{CD<:ConstraintDelegate} = ConstraintRef{Model, CD}

function JuMP.add_constraint(m::Model, pc::PolyConstraint, name::String;
                             domain::AbstractSemialgebraicSet=FullSpace(),
                             basis=MonomialBasis, kwargs...)
    delegate = addpolyconstraint!(m, pc.p, pc.set, domain, basis; kwargs...)
    JuMP.ConstraintRef(m, delegate, JuMP.ScalarShape())
end

getdelegate(c::PolyConstraintRef) = c.index
getslack(c::PolyConstraintRef) = getslack(getdelegate(c))
JuMP.result_dual(c::PolyConstraintRef) = JuMP.result_dual(getdelegate(c))

# Macro
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.EqualTo)
    PolyConstraint(p-s.value, ZeroPoly())
end
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.GreaterThan)
    PolyConstraint(p-s.lower, NonNegPoly())
end
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.LessThan)
    PolyConstraint(s.upper-p, NonNegPoly())
end

function JuMP.build_constraint(_error::Function,
                               p::Union{AbstractPolynomialLike,
                                        AbstractMatrix{<:AbstractPolynomialLike}},
                               s)
    PolyConstraint(p, s)
end
# there is already a method for AbstractMatrix in PSDCone in JuMP so we need a
# more specific here to avoid ambiguity
function JuMP.build_constraint(_error::Function,
                               p::AbstractMatrix{<:AbstractPolynomialLike},
                               s::PSDCone)
    PolyConstraint(p, NonNegPolyMatrix())
end
