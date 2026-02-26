module Objective

import MathOptInterface
const MOI = MathOptInterface

using PolyJuMP

include("to_polynomial.jl")
include("to_reformulated_polynomial.jl")

end # module
