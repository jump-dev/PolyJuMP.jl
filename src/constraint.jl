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
const PolyConstraintRef = ConstraintRef{Model, PolyConstraint}

# Responsible for getting slack and dual values
abstract type ConstraintDelegate end

function JuMP.addconstraint(m::Model, pc::PolyConstraint, name::String; domain::AbstractSemialgebraicSet=FullSpace(), basis=MonomialBasis, kwargs...)
    delegates = getdelegates(m)
    delegate = addpolyconstraint!(m, pc.p, pc.set, domain, basis; kwargs...)
    push!(delegates, delegate)
    JuMP.ConstraintRef(m, delegate)
end

getdelegate(c::PolyConstraintRef) = c.index.delegate
getslack(c::PolyConstraintRef) = getslack(getdelegate(c))
JuMP.resultdual(c::PolyConstraintRef) = JuMP.resultdual(getdelegate(c))

# Macro
function JuMP.buildconstraint(_error::Function, p::AbstractPolynomialLike, s::MOI.EqualTo)
    PolyConstraint(p-s.value, ZeroPoly())
end
function JuMP.buildconstraint(_error::Function, p::AbstractPolynomialLike, s::MOI.GreaterThan)
    PolyConstraint(p-s.lower, NonNegPoly())
end
function JuMP.buildconstraint(_error::Function, p::AbstractPolynomialLike, s::MOI.LessThan)
    PolyConstraint(s.upper-p, NonNegPoly())
end

function JuMP.buildconstraint(_error::Function, np::Union{AbstractPolynomialLike, AbstractMatrix{<:AbstractPolynomialLike}}, s)
    PolyConstraint(p, s)
end
# there is already a method for AbstractMatrix in PSDCone in JuMP so we need a more specific here to avoid ambiguity
function JuMP.buildconstraint(_error::Function, p::AbstractMatrix{<:AbstractPolynomialLike}, s::PSDCone)
    PolyConstraint(p, NonNegPolyMatrix())
end
