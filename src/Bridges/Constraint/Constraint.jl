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

include("zero_polynomial.jl")
include("zero_polynomial_in_algebraic_set.jl")
include("plus_minus.jl")
include("to_polynomial.jl")

end # module
