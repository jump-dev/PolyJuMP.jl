module PolyJuMP

import MutableArithmetics
const MA = MutableArithmetics

# MultivariatePolynomials extension

using MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases
using MultivariateMoments
using SemialgebraicSets

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
