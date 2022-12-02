module TestConstraint

using Test

using MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets

using JuMP
using PolyJuMP

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
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ResultCount()) == 1
    @test MOI.get.(model, MOI.VariablePrimal(), vars) ≈ [1.0, √2/2, √2/2]
    #@test MOI.get(model, MOI.ConstraintDual(), c1) ≈ √2
    #@test MOI.get(model, MOI.ConstraintDual(), c2) ≈ √2/2
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
