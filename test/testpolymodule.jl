module TestPolyModule
using JuMP
using PolyJuMP
using MultivariatePolynomials

struct TestNonNegConstraint
    basis
end
TestNonNegConstraint() = TestNonNegConstraint(MonomialBasis)
struct TestNonNegMatrixConstraint
    basis
end
TestNonNegMatrixConstraint() = TestNonNegMatrixConstraint(MonomialBasis)
struct TestConstraint <: PolyJuMP.ConstraintDelegate
    p
    set
    domain
    kwargs
end
PolyJuMP.addpolyconstraint!(m::JuMP.Model, p, s::Union{TestNonNegConstraint, TestNonNegMatrixConstraint}, domain; kwargs...) = TestConstraint(p, s, domain, kwargs)

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, TestNonNegConstraint)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPolyMatrix, TestNonNegMatrixConstraint)
end

end
