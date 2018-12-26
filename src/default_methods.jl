# Free polynomial
function JuMP.variable_type(model::JuMP.AbstractModel, p::Poly)
    return polytype(model, p, p.polynomial_basis)
end
function polytype(model::JuMP.AbstractModel, ::Poly, pb::AbstractPolynomialBasis)
    return MultivariatePolynomials.polynomialtype(pb, JuMP.VariableRef)
end

function JuMP.add_variable(model::JuMP.AbstractModel, v::Variable{<:Poly},
                           name::String="")
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
    return polynomial(_newvar, v.p.polynomial_basis)
end

# NonNegPoly and NonNegPolyMatrix
function addpolyconstraint!(model::JuMP.Model, p,
                            s::Union{NonNegPoly, NonNegPolyMatrix},
                            domain, basis; kwargs...)
    return addpolyconstraint!(model, p, getdefault(model, s), domain, basis;
                              kwargs...)
end

# ZeroPoly
struct ZeroConstraint{MT <: AbstractMonomial, MVT <: AbstractVector{MT}, F <: MOI.AbstractVectorFunction} <: ConstraintDelegate
    # F is typically VectorAffineFunction or VectorQuadraticFunction
    zero_constraints::JuMP.ConstraintRef{JuMP.Model,MOI.ConstraintIndex{F,MOI.Zeros}}
    x::MVT
end

JuMP.dual(c::ZeroConstraint) = measure(JuMP.dual.(c.zero_constraints), c.x)

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
