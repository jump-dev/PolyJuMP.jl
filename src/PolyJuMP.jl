module PolyJuMP

using MultivariatePolynomials
using MultivariateMoments
using SemialgebraicSets

# MultivariatePolynomials extension

include("basis.jl")

# MOI extension

using MathOptInterface
const MOI = MathOptInterface

include("zero_polynomial.jl")

# JuMP extension

using JuMP

include("variable.jl")
include("constraint.jl")

include("data.jl")
include("default.jl")

end # module
