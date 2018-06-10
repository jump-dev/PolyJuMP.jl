using JuMP
using PolyJuMP
using Base.Test

using SemialgebraicSets

using MultivariatePolynomials
using DynamicPolynomials
#using TypedPolynomials

include("testpolymodule.jl")

include("polymodule.jl")
include("variable.jl")
include("constraint.jl")

if isdir(Pkg.dir("SumOfSquares"))
    #include("sumofsquares.jl")
end
