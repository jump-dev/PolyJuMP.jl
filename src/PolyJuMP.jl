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

# Bridges
const MOIB = MOI.Bridges
include("zero_polynomial_bridge.jl")
include("zero_polynomial_in_algebraic_set_bridge.jl")
include("zero_polynomial_in_semialgebraic_set_bridge.jl")

# JuMP extension

using JuMP

include("variable.jl")
include("constraint.jl")

include("data.jl")
include("default.jl")

end # module
