using MultivariatePolynomials
using JuMP
using PolyJuMP
using Base.Test

module TestPolyModule
using JuMP
type TestPoly
    nonnegative::Bool
    monotype::Symbol
    x
end
createpoly(m::JuMP.Model, monotype::Symbol, x) = TestPoly(false, monotype, x)
createnonnegativepoly(m::JuMP.Model, monotype::Symbol, x) = TestPoly(true, monotype, x)
end

include("polymodule.jl")
include("variable.jl")

if isdir(Pkg.dir("SumOfSquares"))
    include("sumofsquares.jl")
end
