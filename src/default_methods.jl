# Free polynomial
JuMP.variabletype(m::JuMP.Model, p::Poly) = polytype(m, p, p.polynomial_basis)
function polytype(m::JuMP.Model, ::Poly, pb::AbstractPolynomialBasis)
    MultivariatePolynomials.polynomialtype(pb, JuMP.VariableRef)
end

function createpoly(m::JuMP.Model, p::Poly, binary::Bool, integer::Bool)
    function _newvar(i)
        v = VariableRef(m)
        if binary
            JuMP.set_binary(v)
        end
        if integer
            JuMP.set_integer(v)
        end
        v
    end
    polynomial(_newvar, p.polynomial_basis)
end

# NonNegPoly and NonNegPolyMatrix
addpolyconstraint!(m::JuMP.Model, p, s::Union{NonNegPoly, NonNegPolyMatrix}, domain, basis; kwargs...) = addpolyconstraint!(m, p, getdefault(m, s), domain, basis; kwargs...)

# ZeroPoly
struct ZeroConstraint{MT <: AbstractMonomial, MVT <: AbstractVector{MT}, F <: MOI.AbstractVectorFunction} <: ConstraintDelegate
    # F is typically VectorAffineFunction or VectorQuadraticFunction
    zero_constraints::JuMP.ConstraintRef{JuMP.Model,MOI.ConstraintIndex{F,MOI.Zeros}}
    x::MVT
end
if VERSION < v"0.7-"
    function ZeroConstraint(zero_constraints::JuMP.ConstraintRef{JuMP.Model,MOI.ConstraintIndex{F,MOI.Zeros}}, x::MVT) where {MT <: AbstractMonomial, MVT <: AbstractVector{MT}, F <: MOI.AbstractVectorFunction}
        ZeroConstraint{MT, MVT, F}(zero_constraints, x)
    end
end

JuMP.result_dual(c::ZeroConstraint) = measure(JuMP.result_dual.(c.zero_constraints), c.x)

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::FullSpace, basis)
    coeffs = collect(coefficients(p))
    c = JuMP.build_constraint(error, coeffs, MOI.Zeros(length(coeffs)))
    zero_constraints = JuMP.add_constraint(m, c)
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
