module TestConstraint

using Test

using MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets

using JuMP
using PolyJuMP

function _test_solution(model, vars, c1, c2)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ResultCount()) == 1
    @test MOI.get.(model, MOI.VariablePrimal(), vars) ≈ [1.0, √2/2, √2/2]
    @test MOI.get(model, MOI.ConstraintDual(), c1) ≈ √2
    @test MOI.get(model, MOI.ConstraintDual(), c2) ≈ √2/2
end

function test_algebraic(var, T)
    model = PolyJuMP.AlgebraicKKT.Optimizer{T}()
    t = MP.similarvariable(var, Val{:t})
    x = MP.similarvariable(var, Val{:x})
    y = MP.similarvariable(var, Val{:y})
    z = zero(T)
    o = one(T)
    vars = MOI.add_variables(model, 3)
    c1 = MOI.add_constraint(model, PolyJuMP.ScalarPolynomialFunction(o * t + z, vars[1:1]), MOI.LessThan(o))
    c2 = MOI.add_constraint(model, PolyJuMP.ScalarPolynomialFunction(o * x^2 + o * y^2 - o * t^2, vars), MOI.LessThan(z))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj = PolyJuMP.ScalarPolynomialFunction(o * x + o * y, vars[2:3])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    _test_solution(model, vars, c1, c2)
end

function test_quad(var, T)
    model = MOI.instantiate(PolyJuMP.AlgebraicKKT.Optimizer{T}, with_bridge_type=T)
    z = zero(T)
    o = one(T)
    vars = MOI.add_variables(model, 3)
    c1 = MOI.add_constraint(model, o * vars[1] + z, MOI.LessThan(o))
    c2 = MOI.add_constraint(model, o * vars[2] * vars[2] + o * vars[3] * vars[3] - o * vars[1] * vars[1], MOI.LessThan(z))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    t = MP.similarvariable(var, Val{:t})
    x = MP.similarvariable(var, Val{:x})
    y = MP.similarvariable(var, Val{:y})
    obj = PolyJuMP.ScalarPolynomialFunction(o * x + o * y, vars[2:3])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    _test_solution(model, vars, c1, c2)
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

end

import DynamicPolynomials
DynamicPolynomials.@polyvar(x)
TestConstraint.runtests(x, Float64)
import TypedPolynomials
TypedPolynomials.@polyvar(y)
TestConstraint.runtests(y, Float64)
