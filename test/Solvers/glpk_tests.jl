include("../Tests/Tests.jl")

using Test
using JuMP
using GLPK
optimizer = GLPK.Optimizer()
bridged = MOI.Bridges.full_bridge_optimizer(optimizer, Float64)
config = MOI.Test.Config(atol = 1e-5, rtol = 1e-5, query = false)
exclude = ["plus_minus"]
Tests.linear_test(bridged, config, exclude)
