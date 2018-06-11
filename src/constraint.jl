export ZeroPoly, NonNegPoly, NonNegPolyMatrix
export getslack

abstract type PolynomialSet end
struct ZeroPoly{B<:AbstractPolynomialBasis} <: PolynomialSet end
struct NonNegPoly <: PolynomialSet end
struct NonNegPolyMatrix <: PolynomialSet end

struct PolyConstraint{PT, ST<:PolynomialSet} <: JuMP.AbstractConstraint
    p::PT # typically either be a polynomial or a Matrix of polynomials
    set::ST
end
const PolyConstraintRef = ConstraintRef{Model, PolyConstraint}

# Responsible for getting slack and dual values
abstract type ConstraintDelegate end

function JuMP.addconstraint(m::Model, pc::PolyConstraint; domain::AbstractSemialgebraicSet=FullSpace(), kwargs...)
    delegates = getdelegates(m)
    delegate = addpolyconstraint!(m, pc.p, pc.set, domain; kwargs...)
    push!(delegates, delegate)
    m.internalModelLoaded = false
    PolyConstraintRef(m, length(delegates))
end

getdelegate(c::PolyConstraintRef) = getdelegates(c.m)[c.idx]
getslack(c::PolyConstraintRef) = getslack(getdelegate(c))
JuMP.getdual(c::PolyConstraintRef) = getdual(getdelegate(c))

# Macro
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
