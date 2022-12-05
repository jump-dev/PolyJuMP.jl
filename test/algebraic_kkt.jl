module TestAlgebraicKKT

using Test

using MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets

using JuMP
using PolyJuMP

_optimize!(model::MOI.ModelLike) = MOI.optimize!(model)
_optimize!(model::Model) = optimize!(model)
_dual(model::MOI.ModelLike, c) = MOI.get(model, MOI.ConstraintDual(), c)
_dual(model::Model, c) = dual(c)

function _test_solution(model, vars, c1, c2)
    _optimize!(model)
    @test MOI.get(model, MOI.ResultCount()) == 1
    @test MOI.get.(model, MOI.VariablePrimal(), vars) ≈ [1.0, √2 / 2, √2 / 2]
    @test _dual(model, c1) ≈ √2
    @test _dual(model, c2) ≈ √2 / 2
end

function test_algebraic(var, T)
    model = PolyJuMP.AlgebraicKKT.Optimizer{T}()
    t = MP.similarvariable(var, Val{:t})
    x = MP.similarvariable(var, Val{:x})
    y = MP.similarvariable(var, Val{:y})
    z = zero(T)
    o = one(T)
    vars = MOI.add_variables(model, 3)
    c1 = MOI.add_constraint(
        model,
        PolyJuMP.ScalarPolynomialFunction(o * t + z, vars[1:1]),
        MOI.LessThan(o),
    )
    c2 = MOI.add_constraint(
        model,
        PolyJuMP.ScalarPolynomialFunction(o * x^2 + o * y^2 - o * t^2, vars),
        MOI.LessThan(z),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj = PolyJuMP.ScalarPolynomialFunction(o * x + o * y, vars[2:3])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    return _test_solution(model, vars, c1, c2)
end

function _test_linquad(T, F1, O)
    model = MOI.instantiate(
        PolyJuMP.AlgebraicKKT.Optimizer{T},
        with_bridge_type = T,
    )
    z = zero(T)
    o = one(T)
    vars = MOI.add_variables(model, 3)
    if F1 === MOI.VariableIndex
        c1 = MOI.add_constraint(model, vars[1], MOI.LessThan(o))
    elseif F1 == MOI.ScalarAffineFunction
        c1 = MOI.add_constraint(model, o * vars[1] + z, MOI.LessThan(o))
    end
    c2 = MOI.add_constraint(
        model,
        o * vars[2] * vars[2] + o * vars[3] * vars[3] - o * vars[1] * vars[1],
        MOI.LessThan(z),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    if O === MOI.ScalarAffineFunction
        obj = o * vars[2] + o * vars[3]
    elseif O == MOI.ScalarQuadraticFunction
        obj = o * vars[2] + o * vars[3] + z * vars[2] * vars[3]
    end
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    return _test_solution(model, vars, c1, c2)
end

function test_linquad(var, T)
    for F1 in [MOI.VariableIndex, MOI.ScalarAffineFunction]
        for O in [MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction]
            _test_linquad(T, F1, O)
        end
    end
end

function _test_JuMP(F1, O)
    model = Model(PolyJuMP.AlgebraicKKT.Optimizer)
    @variable(model, t)
    @variable(model, x)
    @variable(model, y)
    if F1 === MOI.VariableIndex
        @constraint(model, c1, t in MOI.LessThan(1.0))
    elseif F1 == MOI.ScalarAffineFunction
        @constraint(model, c1, t <= 1)
    elseif F1 == 1
        @NLconstraint(model, c1, t <= 1)
    elseif F1 == 2
        @NLconstraint(model, c1, 1 * t <= 1)
    end
    @constraint(model, c2, x^2 + y^2 <= t^2)
    if O === MOI.ScalarAffineFunction
        @objective(model, Max, x + y)
    elseif O == MOI.ScalarQuadraticFunction
        @objective(model, Max, x + y + 0 * x * y)
    elseif O == 1
        @NLobjective(model, Max, x + y)
    elseif O == 2
        @NLobjective(model, Max, x - (-y))
    elseif O == 3
        @NLobjective(model, Max, -(-x) + y)
    else
        @NLobjective(model, Max, 1 * x + 1 * y)
    end
    return _test_solution(model, [t, x, y], c1, c2)
end

function test_JuMP(var, T)
    for F1 in [MOI.VariableIndex, MOI.ScalarAffineFunction, 1, 2]
        for O in
            [MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction, 1, 2, 3, 4]
            _test_JuMP(F1, O)
        end
    end
end

function runtests(var, T)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)(var, T)
            end
        end
    end
end

end # module

import DynamicPolynomials
DynamicPolynomials.@polyvar(x)
TestAlgebraicKKT.runtests(x, Float64)
import TypedPolynomials
TypedPolynomials.@polyvar(y)
TestAlgebraicKKT.runtests(y, Float64)
