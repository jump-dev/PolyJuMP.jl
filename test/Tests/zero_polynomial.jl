using Test
using JuMP
using MultivariateMoments
using PolyJuMP
using DynamicPolynomials

function zero_polynomial_test(optimizer::MOI.AbstractOptimizer,
                              config::MOI.Test.TestConfig)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)
    @variable(model, β ≤ 1)
    @variable(model, γ)

    @polyvar x y
    cref = @constraint(model, (α - β) * x * y == 0)
    cγ = @constraint(model, γ * x * y in PolyJuMP.ZeroPoly())

    @objective(model, Max, α + γ)
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

    for μ in [dual(cref), moments(cref), dual(cγ), moments(cγ)]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 1
        @test moment_value(moments(μ)[1]) ≈ -1.0 atol=atol rtol=rtol
        @test monomial(moments(μ)[1]) == x*y
    end

    F = MOI.VectorAffineFunction{Float64}
    S = PolyJuMP.ZeroPolynomialSet{FullSpace,MonomialBasis,Monomial{true},
                                   MonomialVector{true}}
    @test Set(MOI.get(model, MOI.ListOfConstraints())) == Set([
        (MOI.SingleVariable, MOI.LessThan{Float64}), (F, S), (MOI.VectorOfVariables, S)])
    @testset "Delete" begin
        test_delete_bridge(model, cref, 3, ((F, MOI.Zeros, 0),))
    end
end

linear_tests["zero_polynomial"] = zero_polynomial_test
