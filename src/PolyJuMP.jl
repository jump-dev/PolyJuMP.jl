module PolyJuMP

import MutableArithmetics as MA

# MultivariatePolynomials extension

import MultivariatePolynomials as MP
import MultivariateBases as MB
import MultivariateMoments as MM
import SemialgebraicSets as SS

# MOI extension

import MathOptInterface as MOI

include("attributes.jl")
include("zero_polynomial.jl")
include("functions.jl")

# Bridges
include("Bridges/Bridges.jl")
include("nl_to_polynomial.jl")

# JuMP extension

using JuMP

include("variable.jl")
include("constraint.jl")

include("data.jl")
include("default.jl")

include("KKT/KKT.jl")

end # module
