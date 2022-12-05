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

using MathOptInterface
const MOI = MathOptInterface

include("attributes.jl")
include("zero_polynomial.jl")
include("functions.jl")

# Bridges
include("Bridges/Bridges.jl")

# JuMP extension

using JuMP

include("variable.jl")
include("constraint.jl")

include("data.jl")
include("default.jl")

include("AlgebraicKKT/AlgebraicKKT.jl")

end # module
