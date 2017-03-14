using MultivariatePolynomials
using JuMP
using PolyJuMP
using Base.Test

module TestPolyModule
using MultivariatePolynomials
using JuMP
type TestPoly
    nonnegative::Bool
    monotype::Symbol
    x
end
type TestCon
    nonnegative::Bool
    p::Polynomial
    domain
end
createpoly(m::JuMP.Model, monotype::Symbol, x) = TestPoly(false, monotype, x)
createnonnegativepoly(m::JuMP.Model, monotype::Symbol, x) = TestPoly(true, monotype, x)
#addpolyeqzeroconstraint(m::JuMP.Model, p, domain::AlgebraicSet) = TestCon(false, p, domain)
#addpolynonnegativeconstraint(m::JuMP.Model, p, domain::BasicSemialgebraicSet) = TestCon(true, p, domain)
end

include("polymodule.jl")
include("variable.jl")
include("constraint.jl")

if isdir(Pkg.dir("SumOfSquares"))
    #include("sumofsquares.jl")
end
