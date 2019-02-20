using Test
using JuMP
using SemialgebraicSets
using PolyJuMP
using DynamicPolynomials

function zero_polynomial_in_algebraic_set_test(optimizer::MOI.AbstractOptimizer,
                                               config::MOI.Test.TestConfig)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    model = JuMP.direct_model(optimizer)

    @variable(model, α)
    @variable(model, β ≤ 1)

    @polyvar x y
    @constraint(model, α * x - β * y == 0, domain = @set x == y)

    @objective(model, Max, α)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 1.0 atol=atol rtol=rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 1.0 atol=atol rtol=rtol
    @test value(β) ≈ 1.0 atol=atol rtol=rtol
    @test value(UpperBoundRef(β)) ≈ 1.0 atol=atol rtol=rtol

    @test dual_status(model) == MOI.FEASIBLE_POINT
    @test dual(UpperBoundRef(β)) ≈ -1.0 atol=atol rtol=rtol
end
