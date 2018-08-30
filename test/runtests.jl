using Compat
using Compat.Test

using MultivariatePolynomials
using SemialgebraicSets

using DynamicPolynomials
#using TypedPolynomials

using JuMP
using PolyJuMP

include("utilities.jl")
include("testpolymodule.jl")

include("polymodule.jl")
include("variable.jl")
include("constraint.jl")

#if isdir(Pkg.dir("SumOfSquares"))
#    #include("sumofsquares.jl")
#end
