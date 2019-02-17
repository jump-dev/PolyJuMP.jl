using Test
using JuMP
using SemialgebraicSets
# `MultivariateMoments.value` clashes with `JuMP.value`
import MultivariateMoments
const MM = MultivariateMoments
using PolyJuMP
using DynamicPolynomials

function plus_minus_test(optimizer::MOI.AbstractOptimizer,
                         config::MOI.Test.TestConfig;
                         polymodule = nothing)
    atol = config.atol
    rtol = config.rtol

    MOI.empty!(optimizer)
    model = JuMP.direct_model(optimizer)

    if polymodule !== nothing
        setpolymodule!(model, polymodule)
    end

    @variable(model, α)
    @variable(model, β ≤ 1)

    @polyvar x y
    # x^2 ≤ 0 implies that x = 0
    #  VectorAffineFunction
    @constraint(model, α * (x + y) - β * y == 0, domain = @set x^2 ≤ 0)
    #  VectorOfVariables
    @constraint(model, α * x in PolyJuMP.ZeroPoly(),
                domain = @set x^2 ≤ 0)

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
