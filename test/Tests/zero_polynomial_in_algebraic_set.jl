if !@isdefined(linear_tests)
    include("utilities.jl")
    const linear_tests = Dict{String,Function}()
end

import MultivariateBases as MB
import SparseArrays
import StarAlgebras as SA
using Test
using JuMP
using SemialgebraicSets
using MultivariateMoments
using PolyJuMP
using DynamicPolynomials

function _zero_polynomial_in_algebraic_set_test(
    optimizer,
    config::MOI.Test.Config,
    x,
    y,
    domain,
    domain_type = typeof(domain),
)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)
    @variable(model, β ≤ 1)

    cref = @constraint(model, α * x - β * y == 0, domain = domain)

    @objective(model, Max, α)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 1.0 atol = atol rtol = rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 1.0 atol = atol rtol = rtol
    @test value(β) ≈ 1.0 atol = atol rtol = rtol
    @test value(UpperBoundRef(β)) ≈ 1.0 atol = atol rtol = rtol

    @test dual_status(model) == MOI.FEASIBLE_POINT
    @test dual(UpperBoundRef(β)) ≈ -1.0 atol = atol rtol = rtol

    μ = dual(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 2
    @test moment_value(moments(μ)[1]) ≈ -1.0 atol = atol rtol = rtol
    @test moments(μ)[1].polynomial == MB.BasisElement{MB.Monomial}(x^0 * y)
    @test moment_value(moments(μ)[2]) ≈ -1.0 atol = atol rtol = rtol
    @test moments(μ)[2].polynomial == MB.BasisElement{MB.Monomial}(x * y^0)

    # The moments are given for the monomials of the polynomial and not only
    # for the standard monomials of the quotient
    μ = moments(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 2
    @test moment_value(moments(μ)[1]) ≈ -1.0 atol = atol rtol = rtol
    @test moments(μ)[1].polynomial == MB.BasisElement{MB.Monomial}(x^0 * y)
    @test moment_value(moments(μ)[2]) ≈ -1.0 atol = atol rtol = rtol
    @test moments(μ)[2].polynomial == MB.BasisElement{MB.Monomial}(x * y^0)

    F = MOI.VectorAffineFunction{Float64}
    full_basis_type = typeof(MB.FullBasis{MB.Monomial}(x))
    sub_basis_type = MB.explicit_basis_type(full_basis_type)
    S = PolyJuMP.ZeroPolynomialSet{domain_type,full_basis_type,sub_basis_type}
    @test MOI.get(model, MOI.ListOfConstraintTypesPresent()) ==
          [(MOI.VariableIndex, MOI.LessThan{Float64}), (F, S)]
    @testset "Delete" begin
        quotient_basis_type = MB.QuotientBasis{
            eltype(sub_basis_type),
            SA.key_type(sub_basis_type),
            sub_basis_type,
            SparseArrays.SparseMatrixCSC{Float64,Int},
        }
        ST = PolyJuMP.ZeroPolynomialSet{
            FullSpace,
            full_basis_type,
            quotient_basis_type,
        }
        test_delete_bridge(model, cref, 2, ((F, MOI.Zeros, 0), (F, ST, 0)))
    end
end

function zero_polynomial_in_algebraic_set_test(
    optimizer,
    config::MOI.Test.Config,
)
    @polyvar x y
    return _zero_polynomial_in_algebraic_set_test(
        optimizer,
        config,
        x,
        y,
        @set(x == y),
    )
end

linear_tests["zero_polynomial_in_algebraic_set"] =
    zero_polynomial_in_algebraic_set_test

# The equation is redundant so that it is not a Gröbner basis. The Gröbner basis
# is not computed and the Macaulay matrix is used instead.
function zero_polynomial_in_algebraic_set_macaulay_test(
    optimizer,
    config::MOI.Test.Config,
)
    @polyvar x y
    domain = @set x == y && 2x == 2y
    # The domain is truncated at the degree of the polynomial and equalities
    _zero_polynomial_in_algebraic_set_test(
        optimizer,
        config,
        x,
        y,
        domain,
        typeof(SemialgebraicSets.truncated(domain, 1)),
    )
    @test !SemialgebraicSets.is_gröbner_basis(ideal(domain))
    return
end

linear_tests["zero_polynomial_in_algebraic_set_macaulay"] =
    zero_polynomial_in_algebraic_set_macaulay_test
