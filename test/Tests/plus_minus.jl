using Test
using JuMP
using SemialgebraicSets
using PolyJuMP
using DynamicPolynomials

function plus_minus_test(optimizer,
                         config::MOI.Test.Config;
                         polymodule = nothing)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    if polymodule !== nothing
        setpolymodule!(model, polymodule)
    end

    @variable(model, α)
    @variable(model, β ≤ 1)

    @polyvar x y
    # x^2 ≤ 0 implies that x = 0
    #  VectorAffineFunction
    c1 = @constraint(model, α * (x + y) - β * y == 0, domain = @set x^2 ≤ 0)
    #  VectorOfVariables
    c2 = @constraint(model, α * x in PolyJuMP.ZeroPoly(),
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

    _NonNegType(::ConstraintRef{<:JuMP.AbstractModel, MOI.ConstraintIndex{MOI.VectorOfVariables, PolyJuMP.PlusMinusSet{S}}}) where S = S

    @testset "Delete" begin
        F = MOI.VectorOfVariables
        G = MOI.VectorAffineFunction{Float64}
        NonNeg = _NonNegType(c2)
        S = PolyJuMP.PlusMinusSet{NonNeg}
        @test Set(MOI.get(model, MOI.ListOfConstraintTypesPresent())) == Set([
            (MOI.VariableIndex, MOI.LessThan{Float64}), (G, S), (F, S)])
        test_delete_bridge(model, c2, 2, ((F, NonNeg, 0), (G, NonNeg, 0)))
        test_delete_bridge(model, c1, 2, ((F, NonNeg, 0), (G, NonNeg, 0)))
    end
end

linear_tests["plus_minus"] = plus_minus_test
