using Test

using MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases
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
include("functions.jl")

include("Mock/mock_tests.jl")
