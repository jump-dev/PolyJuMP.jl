using Test

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

include("zero_polynomial_bridge.jl")
include("zero_polynomial_in_algebraic_set_bridge.jl")
