export ZeroPoly, NonNegPoly, NonNegPolyMatrix
export getslack

abstract type PolynomialSet end
struct ZeroPoly <: PolynomialSet end
struct NonNegPoly <: PolynomialSet end
struct NonNegPolyMatrix <: PolynomialSet end

struct Constraint{PT, ST<:PolynomialSet} <: JuMP.AbstractConstraint
    p::PT # typically either be a polynomial or a Matrix of polynomials
    set::ST
end

# Responsible for getting slack and dual values
abstract type ConstraintDelegate end

const ConstraintRef{CD<:ConstraintDelegate} = JuMP.ConstraintRef{Model, CD}

Base.show(io::IO, cref::ConstraintRef) = print(io, "PolyJuMP constraint")

function JuMP.add_constraint(m::Model, pc::Constraint, name::String="";
                             domain::AbstractSemialgebraicSet=FullSpace(),
                             basis=MonomialBasis, kwargs...)
    delegate = addpolyconstraint!(m, pc.p, pc.set, domain, basis; kwargs...)
    JuMP.ConstraintRef(m, delegate, JuMP.ScalarShape())
end

getdelegate(c::ConstraintRef) = c.index
getslack(c::ConstraintRef) = getslack(getdelegate(c))
JuMP.dual(c::ConstraintRef) = JuMP.dual(getdelegate(c))

# Macro
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.EqualTo)
    Constraint(p-s.value, ZeroPoly())
end
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.GreaterThan)
    Constraint(p-s.lower, NonNegPoly())
end
function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::MOI.LessThan)
    Constraint(s.upper-p, NonNegPoly())
end

function JuMP.build_constraint(_error::Function,
                               p::Union{AbstractPolynomialLike,
                                        AbstractMatrix{<:AbstractPolynomialLike}},
                               s)
    Constraint(p, s)
end
# there is already a method for AbstractMatrix in PSDCone in JuMP so we need a
# more specific here to avoid ambiguity
function JuMP.build_constraint(_error::Function,
                               p::AbstractMatrix{<:AbstractPolynomialLike},
                               s::PSDCone)
    Constraint(p, NonNegPolyMatrix())
end
