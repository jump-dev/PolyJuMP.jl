using JuMP
using PolyJuMP

using Compat
using Compat.Test

using SemialgebraicSets

using MultivariatePolynomials
using DynamicPolynomials
#using TypedPolynomials

include("utilities.jl")
include("testpolymodule.jl")

include("polymodule.jl")
include("variable.jl")
include("constraint.jl")

#if isdir(Pkg.dir("SumOfSquares"))
#    include("sumofsquares.jl")
#end
