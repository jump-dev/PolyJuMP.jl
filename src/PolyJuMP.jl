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
const MOIB = MOI.Bridges
include("zero_polynomial_bridge.jl")
include("zero_polynomial_in_algebraic_set_bridge.jl")
include("plus_minus_bridge.jl")
include("to_polynomial_bridge.jl")

# JuMP extension

using JuMP

include("variable.jl")
include("constraint.jl")

include("data.jl")
include("default.jl")

include("AlgebraicKKT/AlgebraicKKT.jl")

end # module
