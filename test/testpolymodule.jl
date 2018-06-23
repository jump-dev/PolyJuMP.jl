module TestPolyModule
using JuMP
using PolyJuMP
using MultivariatePolynomials

struct TestNonNegConstraint <: PolyJuMP.PolynomialSet end
struct TestNonNegMatrixConstraint <: PolyJuMP.PolynomialSet end
struct TestConstraint <: PolyJuMP.ConstraintDelegate
    p
    set
    domain
    basis
    kwargs
end
PolyJuMP.addpolyconstraint!(m::JuMP.Model, p, s::Union{TestNonNegConstraint, TestNonNegMatrixConstraint}, domain, basis; kwargs...) = TestConstraint(p, s, domain, basis, kwargs)

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, TestNonNegConstraint)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPolyMatrix, TestNonNegMatrixConstraint)
end

end
