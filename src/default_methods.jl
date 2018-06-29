# Free polynomial
JuMP.variabletype(m::JuMP.Model, p::Poly) = polytype(m, p, p.polynomial_basis)
function polytype(m::JuMP.Model, ::Poly, pb::AbstractPolynomialBasis)
    MultivariatePolynomials.polynomialtype(pb, JuMP.VariableRef)
end

function createpoly(m::JuMP.Model, p::Poly, binary::Bool, integer::Bool)
    function _newvar(i)
        v = VariableRef(m)
        if binary
            JuMP.setbinary(v)
        end
        if integer
            JuMP.setinteger(v)
        end
        v
    end
    polynomial(_newvar, p.polynomial_basis)
end

# NonNegPoly and NonNegPolyMatrix
addpolyconstraint!(m::JuMP.Model, p, s::Union{NonNegPoly, NonNegPolyMatrix}, domain, basis; kwargs...) = addpolyconstraint!(m, p, getdefault(m, s), domain, basis; kwargs...)

# ZeroPoly
struct ZeroConstraint{MT <: AbstractMonomial, MVT <: AbstractVector{MT}, JC <: JuMP.AbstractConstraint} <: ConstraintDelegate
    zero_constraints::Vector{JuMP.ConstraintRef{JuMP.Model, JC}} # These are typically affine or quadratic equality constraints
    x::MVT
end
function ZeroConstraint(zero_constraints::Vector{JuMP.ConstraintRef{JuMP.Model, JC}}, x::MVT) where {MT <: AbstractMonomial, MVT <: AbstractVector{MT}, JC <: JuMP.AbstractConstraint}
    ZeroConstraint{MT, MVT, JC}(zero_constraints, x)
end

JuMP.resultdual(c::ZeroConstraint) = measure(JuMP.resultdual.(c.zero_constraints), c.x)

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::FullSpace, basis)
    constraints = JuMP.constructconstraint!.(coefficients(p, basis), :(==))
    zero_constraints = JuMP.addVectorizedConstraint(m, constraints)
    ZeroConstraint(zero_constraints, monomials(p))
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::AbstractAlgebraicSet, basis)
    addpolyconstraint!(m, rem(p, ideal(domain)), s, FullSpace(), basis)
end

struct ZeroConstraintWithDomain{DT<:ConstraintDelegate} <: ConstraintDelegate
    lower::DT
    upper::DT
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::BasicSemialgebraicSet, basis)
    lower = addpolyconstraint!(m,  p, NonNegPoly(), domain, basis)
    upper = addpolyconstraint!(m, -p, NonNegPoly(), domain, basis)
    ZeroConstraintWithDomain(lower, upper)
end
