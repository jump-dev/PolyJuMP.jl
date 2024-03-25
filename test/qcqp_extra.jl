module TestQCQPExtra

using Test

import MathOptInterface as MOI
import MultivariatePolynomials as MP
import PolyJuMP
import JuMP
import Random

MOI.Utilities.@model(
    Model,
    (),
    (MOI.LessThan, MOI.GreaterThan, MOI.EqualTo, MOI.Interval),
    (),
    (),
    (),
    (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
    (),
    (),
)

function MOI.supports(
    ::Model,
    ::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction},
)
    return false
end

function run_tests_e2e()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_e2e")
            @testset "$(name) $T" for T in [Int, Float64]
                getfield(@__MODULE__, name)(T)
            end
        end
    end
    return
end

function run_test_scalar_polynomial_function(xs, samples)
    @testset "qcqp extra $T" for T in [Float64, BigFloat]
        for i in 1:samples
            test_scalar_polynomial_function(xs, T)
        end
    end
    return
end

function run_tests_subs(xs, ys, samples, TS)
    Random.seed!(2024)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_subs")
            @testset "$(name) $T" for T in TS
                for i in 1:samples
                    getfield(@__MODULE__, name)(xs, ys, T)
                end
            end
        end
    end
    return
end

function test_unconstrained_before_projection(T)
    inner = Model{T}()
    optimizer = MOI.Utilities.MockOptimizer(inner)
    model = JuMP.GenericModel{T}(
        () -> PolyJuMP.QCQP.Optimizer{T}(optimizer),
    )
    PolyJuMP.@variable(model, -1 <= a[1:2] <= 1)
    PolyJuMP.@objective(model, Min, a[1]^2 * a[2]^2)
    PolyJuMP.optimize!(model)
    vis = MOI.get(inner, MOI.ListOfVariableIndices())
    @test length(vis) == 4
    F = MOI.ScalarQuadraticFunction{T}
    S = MOI.EqualTo{T}
    cis = MOI.get(inner, MOI.ListOfConstraintIndices{F,S}())
    @test length(cis) == 2
    return
end

function test_unconstrained_after_projection(T)
    inner = Model{T}()
    optimizer = MOI.Utilities.MockOptimizer(inner)
    model = JuMP.GenericModel{T}(
        () -> PolyJuMP.QCQP.Optimizer{T}(optimizer),
    )
    PolyJuMP.@variable(model, -1 <= a <= 1)
    PolyJuMP.@objective(model, Min, a^2)
    PolyJuMP.optimize!(model)
    vis = MOI.get(inner, MOI.ListOfVariableIndices())
    @test length(vis) == 1
    F = MOI.ScalarQuadraticFunction{T}
    S = MOI.EqualTo{T}
    cis = MOI.get(inner, MOI.ListOfConstraintIndices{F,S}())
    @test length(cis) == 0
    return
end

function _random_polynomial(vars, T)
    ms = Random.shuffle(MP.monomials(vars, 1:length(vars)))
    return sum(ms[i] * T(randn()) for i in eachindex(ms) if rand(Bool))
end

function test_subs!_preserves_moi_sync(xs, ys, T)
    p = _random_polynomial(xs, T)
    mois = MOI.VariableIndex.(eachindex(xs))
    vals = T.(randn(length(mois)))
    mask = rand(Bool, length(xs))
    is = Random.shuffle(eachindex(xs)[mask])
    index = Dict{eltype(mois),eltype(xs)}(zip(mois[is], ys[is]))
    moi_map = Dict(zip(xs, mois))
    moivars = [moi_map[v] for v in MP.variables(p)]
    before = PolyJuMP.ScalarPolynomialFunction(p, moivars)
    after, _ = PolyJuMP.QCQP._subs!(before, index)
    bmap = [vals[v.value] for v in before.variables]
    amap = [vals[v.value] for v in after.variables]
    bvalue = before.polynomial(MP.variables(before.polynomial) => bmap)
    avalue = after.polynomial(MP.variables(after.polynomial) => amap)
    # avoid verbose fails
    @test isapprox(Float64(bvalue), Float64(avalue))
    return
end

function test_scalar_polynomial_function(xs, T)
    pick = rand(eachindex(xs))
    ids = Random.shuffle(eachindex(xs))
    poly = sum(T(randn()) * xs[i] for i in ids if i != pick)
    mois = MOI.VariableIndex.(eachindex(xs))
    moi_to_vars = Dict(zip(mois, xs))
    spf = PolyJuMP._scalar_polynomial(moi_to_vars, Any, poly)
    expected = MP.variables(poly)
    actual = [moi_to_vars[m] for m in spf.variables]
    @test length(MP.variables(spf.polynomial)) == length(spf.variables)
    @test expected == actual
    return
end

end # module

using Test

@testset "TestQCQPFinalTouch" begin
    TestQCQPExtra.run_tests_e2e()
end

import DynamicPolynomials
@testset "DynamicPolynomials" begin
    ids = 1:4
    DynamicPolynomials.@polyvar(x[ids])
    DynamicPolynomials.@polyvar(y[ids])
    samples = 10
    types = [Float64, BigFloat] # Rational fails with DynamicPolynomials
    TestQCQPExtra.run_tests_subs(x, y, samples, types)
    TestQCQPExtra.run_test_scalar_polynomial_function(x, samples)
end

import TypedPolynomials
@testset "TypedPolynomials" begin
    TypedPolynomials.@polyvar(z[1:4])
    TypedPolynomials.@polyvar(w[1:4])
    types = [Float64, BigFloat, Rational{BigInt}]
    samples = 10
    TestQCQPExtra.run_tests_subs(z, w, samples, types)
end
