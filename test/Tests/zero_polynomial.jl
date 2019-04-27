using Test
using JuMP
using MultivariateMoments
using PolyJuMP
using DynamicPolynomials

function zero_polynomial_test(optimizer::MOI.AbstractOptimizer,
                              config::MOI.Test.TestConfig)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    model = JuMP.direct_model(optimizer)

    @variable(model, α)
    @variable(model, β ≤ 1)

    @polyvar x y
    cref = @constraint(model, (α - β) * x * y == 0)

    @objective(model, Max, α)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 1.0 atol=atol rtol=rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 1.0 atol=atol rtol=rtol
    @test value(β) ≈ 1.0 atol=atol rtol=rtol
    @test value(UpperBoundRef(β)) ≈ 1.0 atol=atol rtol=rtol
    @test value(cref) isa MultivariatePolynomials.AbstractPolynomial{Float64}
    @test value(cref) ≈ 0.0 * x * y atol=atol rtol=rtol

    @test dual_status(model) == MOI.FEASIBLE_POINT
    @test dual(UpperBoundRef(β)) ≈ -1.0 atol=atol rtol=rtol
    μ = dual(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 1
    @test moment_value(moments(μ)[1]) ≈ -1.0 atol=atol rtol=rtol
    @test monomial(moments(μ)[1]) == x*y
end
