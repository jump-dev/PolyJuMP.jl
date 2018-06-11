__precompile__()

module PolyJuMP

using MultivariatePolynomials
using MultivariateMoments
using SemialgebraicSets
using JuMP

include("basis.jl")

include("variable.jl")
include("constraint.jl")
include("default_methods.jl")

include("data.jl")
include("default.jl")

end # module
