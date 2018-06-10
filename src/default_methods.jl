# Free polynomial
JuMP.variabletype(m::JuMP.Model, p::Poly) = polytype(m, p, p.polynomial_basis)
function polytype(m::JuMP.Model, ::Poly, pb::AbstractPolynomialBasis)
    MultivariatePolynomials.polynomialtype(pb, JuMP.Variable)
end

function createpoly(m::JuMP.Model, p::Poly, category::Symbol)
    polynomial(i -> Variable(m, -Inf, Inf, category), p.polynomial_basis)
end

# NonNegPoly and NonNegPolyMatrix
addpolyconstraint!(m::JuMP.Model, p, s::Union{NonNegPoly, NonNegPolyMatrix}, domain; kwargs...) = addpolyconstraint!(m, p, getdefault(m, s), domain; kwargs...)

# ZeroPoly
struct ZeroConstraint{MT <: AbstractMonomial, MVT <: AbstractVector{MT}, JC <: JuMP.AbstractConstraint} <: ConstraintDelegate
    zero_constraints::Vector{JuMP.ConstraintRef{JuMP.Model, JC}} # These are typically affine or quadratic equality constraints
    x::MVT
end
function ZeroConstraint(zero_constraints::Vector{JuMP.ConstraintRef{JuMP.Model, JC}}, x::MVT) where {MT <: AbstractMonomial, MVT <: AbstractVector{MT}, JC <: JuMP.AbstractConstraint}
    ZeroConstraint{MT, MVT, JC}(zero_constraints, x)
end

JuMP.getdual(c::ZeroConstraint) = measure(getdual.(c.zero_constraints), c.x)

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::FullSpace)
    constraints = JuMP.constructconstraint!.(coefficients(p), :(==))
    zero_constraints = JuMP.addVectorizedConstraint(m, constraints)
    ZeroConstraint(zero_constraints, monomials(p))
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::AbstractAlgebraicSet)
    addpolyconstraint!(m, rem(p, ideal(domain)), s, FullSpace())
end

struct ZeroConstraintWithDomain{DT<:ConstraintDelegate} <: ConstraintDelegate
    lower::DT
    upper::DT
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::BasicSemialgebraicSet)
    lower = addpolyconstraint!(m,  p, NonNegPoly(), domain)
    upper = addpolyconstraint!(m, -p, NonNegPoly(), domain)
    ZeroConstraintWithDomain(lower, upper)
end
