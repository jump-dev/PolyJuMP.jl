using Test
using JuMP
using MultivariateMoments
using PolyJuMP
using DynamicPolynomials
import MultivariatePolynomials as MP

function _zero_polynomial_test(
    optimizer::MOI.AbstractOptimizer,
    config::MOI.Test.Config,
    complex::Bool,
)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)
    @variable(model, β ≤ 1)
    @variable(model, γ)

    @polyvar x y
    z = complex ? 0 + 0im : 0
    cref = @constraint(model, (α - β) * x * y == z)
    p = γ * x * y
    cγ = @constraint(model, p in PolyJuMP.ZeroPoly())

    @objective(model, Max, α + γ)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 1.0 atol = atol rtol = rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 1.0 atol = atol rtol = rtol
    @test value(β) ≈ 1.0 atol = atol rtol = rtol
    @test value(γ) ≈ 0.0 atol = atol rtol = rtol
    @test value(UpperBoundRef(β)) ≈ 1.0 atol = atol rtol = rtol
    @test value(cref) isa MP.AbstractPolynomial{Float64}
    @test value(cref) ≈ 0.0 * x * y atol = atol rtol = rtol
    @test value(cγ) ≈ 0.0 * x * y atol = atol rtol = rtol

    @test dual_status(model) == MOI.FEASIBLE_POINT
    @test dual(UpperBoundRef(β)) ≈ -1.0 atol = atol rtol = rtol

    for μ in [dual(cref), moments(cref), dual(cγ), moments(cγ)]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 1
        @test moment_value(moments(μ)[1]) ≈ -1.0 atol = atol rtol = rtol
        @test monomial(moments(μ)[1]) == x * y
    end

    F = MOI.VectorAffineFunction{Float64}
    ST = PolyJuMP.ZeroPolynomialSet{
        FullSpace,
        MB.MonomialBasis,
        monomial_type(x),
        monomial_vector_type(x),
    }
    SP = PolyJuMP.ZeroPolynomialSet{
        FullSpace,
        MB.MonomialBasis,
        monomial_type(x),
        monomial_vector_type(x),
    }
    @test Set(MOI.get(model, MOI.ListOfConstraintTypesPresent())) == Set([
        (MOI.VariableIndex, MOI.LessThan{Float64}),
        (F, SP),
        (MOI.VectorOfVariables, ST),
    ])
    @testset "Delete" begin
        test_delete_bridge(model, cref, 3, ((F, MOI.Zeros, 0),))
    end
end

function real_zero_polynomial_test(
    optimizer::MOI.AbstractOptimizer,
    config::MOI.Test.Config,
)
    return _zero_polynomial_test(optimizer, config, false)
end
function complex_zero_polynomial_test(
    optimizer::MOI.AbstractOptimizer,
    config::MOI.Test.Config,
)
    return _zero_polynomial_test(optimizer, config, false)
end

linear_tests["zero_polynomial"] = real_zero_polynomial_test
linear_tests["zero_polynomial"] = complex_zero_polynomial_test
