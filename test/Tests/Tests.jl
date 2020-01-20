module Tests

include("utilities.jl")

const linear_tests = Dict{String, Function}()

import MultivariateBases
const MB = MultivariateBases

include("zero_polynomial.jl")
include("zero_polynomial_in_fixed_variables_set.jl")
include("zero_polynomial_in_algebraic_set.jl")
include("plus_minus.jl")

@test_suite linear

end
