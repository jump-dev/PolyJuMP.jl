using MultivariatePolynomials
using JuMP
using PolyJuMP
using Base.Test

module TestPolyModule
using MultivariatePolynomials
using JuMP
using PolyJuMP
type TestPoly{P}
    monotype::Symbol
    x
    category::Symbol
end
type TestCon
    nonnegative::Bool
    p::Polynomial
    domain
end
polytype{MT}(m::JuMP.Model, p::Poly{false, MT}) = TestPoly{false}
nonnegativepolytype{MT}(m::JuMP.Model, p::Poly{true, MT}) = TestPoly{true}
createpoly{MT}(m::JuMP.Model, p::Poly{false, MT}, category::Symbol) = TestPoly{false}(MT, p.x, category)
createnonnegativepoly{MT}(m::JuMP.Model, p::Poly{true, MT}, category::Symbol) = TestPoly{true}(MT, p.x, category)
#addpolyeqzeroconstraint(m::JuMP.Model, p, domain::AlgebraicSet) = TestCon(false, p, domain)
#addpolynonnegativeconstraint(m::JuMP.Model, p, domain::BasicSemialgebraicSet) = TestCon(true, p, domain)
end

include("polymodule.jl")
include("variable.jl")
include("constraint.jl")

if isdir(Pkg.dir("SumOfSquares"))
    #include("sumofsquares.jl")
end
