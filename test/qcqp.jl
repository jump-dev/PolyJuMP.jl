module TestQCQP

using Test

import MathOptInterface as MOI
import MultivariatePolynomials as MP
import PolyJuMP

function _test_decompose(monos, exps)
    vars = MP.variables(monos)
    M = eltype(monos)
    expected = PolyJuMP.QCQP.DataStructures.OrderedDict{M,M}(
        var => var for var in vars
    )
    for exp in exps
        expected[exp[1]] = exp[2]
    end
    quad = PolyJuMP.QCQP.decompose(monos)
    @test quad == expected
end

function test_decompose(x, y, _)
    _test_decompose([x * y], [x * y => y])
    return _test_decompose(
        [x^2, y^3, x^2 * y^3],
        [x^2 => x, y^2 => y, x^2 * y^3 => y^3, y^3 => y^2],
    )
end

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

function _test_objective_or_constraint(x, y, T, obj::Bool)
    inner = Model{T}()
    optimizer = MOI.Utilities.MockOptimizer(inner)
    model = PolyJuMP.JuMP.GenericModel{T}(() -> PolyJuMP.QCQP.Optimizer{T}(optimizer))
    PolyJuMP.@variable(model, 1 <= a <= 2)
    PolyJuMP.@variable(model, -5 <= b <= 3)
    PolyJuMP.@constraint(model, a + b >= 1)
    if obj
        PolyJuMP.@objective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
    else
        PolyJuMP.@constraint(model, 0 <= a^3 - a^2 + 2a*b - b^2 + b^3 <= 1)
    end
    PolyJuMP.optimize!(model)
    vis = MOI.get(inner, MOI.ListOfVariableIndices())
    @test length(vis) == 4
    a, b, b2, a2 = vis
    @test MOI.Utilities.get_bounds(inner, T, a) == (1, 2)
    @test MOI.Utilities.get_bounds(inner, T, b) == (-5, 3)
    @test MOI.Utilities.get_bounds(inner, T, a2) == (1, 4)
    @test MOI.Utilities.get_bounds(inner, T, b2) == (0, 25)
    F = MOI.ScalarQuadraticFunction{T}
    S = MOI.EqualTo{T}
    cis = MOI.get(inner, MOI.ListOfConstraintIndices{F,S}())
    @test length(cis) == 2
    o = one(T)
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[1]) ≈ b2 - o * b * b
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[2]) ≈ a2 - o * a * a
    for ci in cis
        @test MOI.get(inner, MOI.ConstraintSet(), ci) == MOI.EqualTo(zero(T))
    end
    exp = -o * a2 + 2o * a * b - o * b2 + o * a * a2 + o * b * b2
    if obj
        @test MOI.get(inner, MOI.ObjectiveFunction{F}()) ≈ exp
    else
        S = MOI.Interval{T}
        cis = MOI.get(inner, MOI.ListOfConstraintIndices{F,S}())
        @test length(cis) == 1
        @test MOI.get(inner, MOI.ConstraintFunction(), cis[1]) ≈ exp
    end
end

function test_objective(x, y, T)
    return _test_objective_or_constraint(x, y, T, true)
end

function test_constraint(x, y, T)
    return _test_objective_or_constraint(x, y, T, false)
end

function test_objective_and_constraint(x, y, T)
    inner = Model{T}()
    optimizer = MOI.Utilities.MockOptimizer(inner)
    model = PolyJuMP.JuMP.GenericModel{T}(() -> PolyJuMP.QCQP.Optimizer{T}(optimizer))
    PolyJuMP.@variable(model, -2 <= a <= 3)
    PolyJuMP.@variable(model, 5 <= b <= 7)
    PolyJuMP.@constraint(model, 0 <= a^3 <= 1)
    PolyJuMP.@constraint(model, 0 <= b^3 <= 1)
    PolyJuMP.@constraint(model, 0 <= a^3 * b^3 + a^6 <= 1)
    PolyJuMP.@objective(model, Max, a^6 * b^3)
    PolyJuMP.optimize!(model)
    vis = MOI.get(inner, MOI.ListOfVariableIndices())
    @test length(vis) == 7
    a, b, b2, a2, a3, b3, a6 = vis
    @test MOI.Utilities.get_bounds(inner, T, a) == (-2, 3)
    @test MOI.Utilities.get_bounds(inner, T, b) == (5, 7)
    @test MOI.Utilities.get_bounds(inner, T, b2) == (25, 49)
    @test MOI.Utilities.get_bounds(inner, T, a2) == (0, 9)
    @test MOI.Utilities.get_bounds(inner, T, a3) == (-18, 27)
    @test MOI.Utilities.get_bounds(inner, T, b3) == (125, 343)
    @test MOI.Utilities.get_bounds(inner, T, a6) == (0, 729)
    F = MOI.ScalarQuadraticFunction{T}
    S = MOI.EqualTo{T}
    cis = MOI.get(inner, MOI.ListOfConstraintIndices{F,S}())
    @test length(cis) == 5
    o = one(T)
    z = zero(T)
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[1]) ≈ b2 - o * b * b
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[2]) ≈ a2 - o * a * a
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[3]) ≈ a3 - o * a2 * a
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[4]) ≈ b3 - o * b2 * b
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[5]) ≈ a6 - o * a3 * a3
    for ci in cis
        @test MOI.get(inner, MOI.ConstraintSet(), ci) == MOI.EqualTo(zero(T))
    end
    @test MOI.get(inner, MOI.ObjectiveFunction{F}()) ≈ o * a6 * b3
    S = MOI.Interval{T}
    cis = MOI.get(inner, MOI.ListOfConstraintIndices{F,S}())
    @test length(cis) == 3
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[1]) ≈ o * a3 + z * a * b
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[2]) ≈ o * b3 + z * a * b
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[3]) ≈ o * a3 * b3 + o * a6
end

function runtests(x, y)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name) $T" for T in [Int, Float64]
                getfield(@__MODULE__, name)(x, y, T)
            end
        end
    end
end

end # module

using Test

import DynamicPolynomials
@testset "DynamicPolynomials" begin
    DynamicPolynomials.@polyvar(x, y)
    TestQCQP.runtests(x, y)
end
import TypedPolynomials
@testset "TypedPolynomials" begin
    TypedPolynomials.@polyvar(x, y)
    TestQCQP.runtests(x, y)
end