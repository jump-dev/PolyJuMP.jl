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
    (MOI.LessThan, MOI.GreaterThan, MOI.EqualTo),
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

function test_objective(x, y, T)
    inner = Model{T}()
    optimizer = MOI.Utilities.MockOptimizer(inner)
    model = PolyJuMP.JuMP.GenericModel{T}(() -> PolyJuMP.QCQP.Optimizer(optimizer))
    PolyJuMP.@variable(model, 0 <= a)
    PolyJuMP.@variable(model, 0 <= b)
    PolyJuMP.@constraint(model, a + b >= 1)
    PolyJuMP.@objective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
    PolyJuMP.optimize!(model)
    vis = MOI.get(inner, MOI.ListOfVariableIndices())
    @test length(vis) == 4
    F = MOI.ScalarQuadraticFunction{T}
    S = MOI.EqualTo{T}
    cis = MOI.get(inner, MOI.ListOfConstraintIndices{F,S}())
    @test length(cis) == 2
end

function runtests(x, y, T)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
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
    TestQCQP.runtests(x, y, Float64)
end
import TypedPolynomials
@testset "TypedPolynomials" begin
    TypedPolynomials.@polyvar(x, y)
    TestQCQP.runtests(x, y, Float64)
end
