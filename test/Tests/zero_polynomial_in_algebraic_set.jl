using Test
using JuMP
using SemialgebraicSets
using MultivariateMoments
using PolyJuMP
using DynamicPolynomials

function zero_polynomial_in_algebraic_set_test(optimizer,
                                               config::MOI.Test.TestConfig)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)
    @variable(model, β ≤ 1)

    @polyvar x y
    cref = @constraint(model, α * x - β * y == 0, domain = @set x == y)

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
    μ = dual(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 2
    @test moment_value(moments(μ)[1]) ≈ -1.0 atol=atol rtol=rtol
    @test monomial(moments(μ)[1]) == x
    @test moment_value(moments(μ)[2]) ≈ -1.0 atol=atol rtol=rtol
    @test monomial(moments(μ)[2]) == y

    F = MOI.VectorAffineFunction{Float64}
    S = PolyJuMP.ZeroPolynomialSet{typeof(@set x == 1),MonomialBasis,
                                   monomialtype(x), monovectype(x)}
    @test MOI.get(model, MOI.ListOfConstraints()) == [
        (MOI.SingleVariable, MOI.LessThan{Float64}), (F, S)]
    @testset "Delete" begin
        ST = PolyJuMP.ZeroPolynomialSet{FullSpace,MonomialBasis,Monomial{true},
                                        MonomialVector{true}}
        test_delete_bridge(model, cref, 2, ((F, MOI.Zeros, 0), (F, ST, 0)))
    end
end

linear_tests["zero_polynomial_in_algebraic_set"] = zero_polynomial_in_algebraic_set_test
