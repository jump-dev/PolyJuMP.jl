module PolyJuMP

import MutableArithmetics as MA
import StarAlgebras as SA

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

include("model.jl")
include("KKT/KKT.jl")
include("QCQP/QCQP.jl")
include("SAGE/SAGE.jl")

end # module
