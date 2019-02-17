include("zero_polynomial_test.jl")
include("zero_polynomial_in_algebraic_set_test.jl")

using GLPK
const optimizer = GLPK.Optimizer()
const bridged = MOI.Bridges.full_bridge_optimizer(optimizer, Float64)
const config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Zero polynomial" begin
    zero_polynomial_test(bridged, config)
end
@testset "Zero polynomial in algebraic set" begin
    zero_polynomial_in_algebraic_set_test(bridged, config)
end
