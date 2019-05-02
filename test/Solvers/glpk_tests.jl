include("../Tests/Tests.jl")

using Test
using MathOptInterface
const MOI = MathOptInterface
using GLPK
const optimizer = GLPK.Optimizer()
const bridged = MOI.Bridges.full_bridge_optimizer(optimizer, Float64)
const config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Zero polynomial" begin
    Tests.zero_polynomial_test(bridged, config)
end
@testset "Zero polynomial in algebraic set" begin
    Tests.zero_polynomial_in_algebraic_set_test(bridged, config)
end
