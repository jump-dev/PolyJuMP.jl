module TestRelativeEntropy

using Test

import MultivariatePolynomials as MP
using SemialgebraicSets

import DynamicPolynomials
import ECOS

using JuMP
using PolyJuMP

function test_motzkin(x, y, T, solver)
    model = Model(solver)
    PolyJuMP.setpolymodule!(model, PolyJuMP.RelativeEntropy)
    motzkin = x^4 * y^2 + x^2 * y^4 + one(T) - 3x^2 * y^2
    @constraint(model, motzkin in PolyJuMP.RelativeEntropy.AGESet(x^2 * y^2))
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test primal_status(model) == MOI.FEASIBLE_POINT
end

import ECOS
const SOLVERS =
    [optimizer_with_attributes(ECOS.Optimizer, MOI.Silent() => true)]

function runtests(x, y, T)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name) $solver)" for solver in SOLVERS
                getfield(@__MODULE__, name)(x, y, T, solver)
            end
        end
    end
end

end # module

using Test

import DynamicPolynomials
@testset "DynamicPolynomials" begin
    DynamicPolynomials.@polyvar(x, y)
    TestRelativeEntropy.runtests(x, y, Float64)
end

import TypedPolynomials
@testset "DynamicPolynomials" begin
    TypedPolynomials.@polyvar(x, y)
    TestRelativeEntropy.runtests(x, y, Float64)
end
