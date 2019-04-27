module Tests

include("utilities.jl")

const linear_tests = Dict{String, Function}()

include("zero_polynomial.jl")
include("zero_polynomial_in_algebraic_set.jl")
include("plus_minus.jl")

@test_suite linear

end
