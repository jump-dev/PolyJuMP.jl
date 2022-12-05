module Constraint

using MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases
import SemialgebraicSets
const SS = SemialgebraicSets
import MultivariateMoments
const MM = MultivariateMoments
import MathOptInterface
const MOI = MathOptInterface

using PolyJuMP

include("zero_polynomial_bridge.jl")
include("zero_polynomial_in_algebraic_set_bridge.jl")
include("plus_minus_bridge.jl")
include("to_polynomial_bridge.jl")

end # module
