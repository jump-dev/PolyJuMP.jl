module PolyJuMP

using MultivariatePolynomials
using MultivariateMoments
using SemialgebraicSets

using MathOptInterface
const MOI = MathOptInterface

using JuMP

include("basis.jl")

include("variable.jl")
include("constraint.jl")
include("default_methods.jl")

include("data.jl")
include("default.jl")

end # module
