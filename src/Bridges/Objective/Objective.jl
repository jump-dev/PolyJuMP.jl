module Objective

import MathOptInterface as MOI

using PolyJuMP

include("to_polynomial.jl")
include("to_reformulated_polynomial.jl")

end # module
