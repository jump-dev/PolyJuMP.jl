module TestRelativeEntropy

using Test

import MultivariatePolynomials as MP
using SemialgebraicSets

import DynamicPolynomials
import ECOS

using JuMP
using PolyJuMP

function _test_motzkin(x, y, T, solver, set, feasible, square, neg)
    model = Model(solver)
    a = square ? x^2 : x
    b = square ? y^2 : y
    PolyJuMP.setpolymodule!(model, PolyJuMP.RelativeEntropy)
    if neg
        motzkin = -a^2 * b - a * b^2 + one(T) - 3a * b
    else
        motzkin = a^2 * b + a * b^2 + one(T) - 3a * b
    end
    con_ref = @constraint(model, motzkin in set)
    optimize!(model)
    if feasible
        @test termination_status(model) == MOI.OPTIMAL
        @test primal_status(model) == MOI.FEASIBLE_POINT
        if set isa Union{
            PolyJuMP.RelativeEntropy.SignomialSAGESet,
            PolyJuMP.RelativeEntropy.PolynomialSAGESet,
        }
            d = PolyJuMP.RelativeEntropy.decomposition(con_ref; tol = 1e-6)
            p = MP.polynomial(d)
            if set isa PolyJuMP.RelativeEntropy.SignomialSAGESet
                @test p ≈ motzkin atol = 1e-6
            else
                for m in MP.monomials(p - motzkin)
                    @test MP.coefficient(p, m) ≈ MP.coefficient(motzkin, m) atol =
                        1e-6
                end
            end
        end
    else
        @test termination_status(model) == MOI.INFEASIBLE
    end
end

function test_motzkin(x, y, T, solver)
    set = PolyJuMP.RelativeEntropy.SignomialAGESet(x^2 * y^2)
    _test_motzkin(x, y, T, solver, set, true, true, false)
    set = PolyJuMP.RelativeEntropy.SignomialAGESet(x * y)
    _test_motzkin(x, y, T, solver, set, true, false, false)
    set = PolyJuMP.RelativeEntropy.SignomialAGESet(x^4 * y^2)
    _test_motzkin(x, y, T, solver, set, false, true, false)
    set = PolyJuMP.RelativeEntropy.SignomialAGESet(x^2 * y)
    _test_motzkin(x, y, T, solver, set, false, false, false)
    set = PolyJuMP.RelativeEntropy.SignomialSAGESet()
    _test_motzkin(x, y, T, solver, set, true, true, false)
    _test_motzkin(x, y, T, solver, set, true, false, false)
    _test_motzkin(x, y, T, solver, set, false, true, true)
    _test_motzkin(x, y, T, solver, set, false, false, true)
    set = PolyJuMP.RelativeEntropy.PolynomialAGESet(x^2 * y^2)
    _test_motzkin(x, y, T, solver, set, true, true, false)
    set = PolyJuMP.RelativeEntropy.PolynomialAGESet(x * y)
    _test_motzkin(x, y, T, solver, set, false, false, false)
    set = PolyJuMP.RelativeEntropy.PolynomialAGESet(x^4 * y^2)
    _test_motzkin(x, y, T, solver, set, false, true, false)
    set = PolyJuMP.RelativeEntropy.PolynomialAGESet(x^2 * y)
    _test_motzkin(x, y, T, solver, set, false, false, false)
    set = PolyJuMP.RelativeEntropy.PolynomialSAGESet()
    _test_motzkin(x, y, T, solver, set, true, true, false)
    set = PolyJuMP.RelativeEntropy.PolynomialSAGESet()
    _test_motzkin(x, y, T, solver, set, false, false, false)
    return
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
