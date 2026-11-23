module Bridges

import MathOptInterface as MOI
using PolyJuMP

include("Utilities.jl")
# Constraint bridges
include("Constraint/Constraint.jl")
# Objective bridges
include("Objective/Objective.jl")

end # module
