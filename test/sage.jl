module TestSAGE

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
    PolyJuMP.setpolymodule!(model, PolyJuMP.SAGE)
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
            PolyJuMP.SAGE.Signomials{Nothing},
            PolyJuMP.SAGE.Polynomials{Nothing},
        }
            d = PolyJuMP.SAGE.decomposition(con_ref; tol = 1e-6)
            p = MP.polynomial(d)
            if set isa PolyJuMP.SAGE.Signomials
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
    set = PolyJuMP.SAGE.Signomials(x^2 * y^2)
    _test_motzkin(x, y, T, solver, set, true, true, false)
    set = PolyJuMP.SAGE.Signomials(x * y)
    _test_motzkin(x, y, T, solver, set, true, false, false)
    set = PolyJuMP.SAGE.Signomials(x^4 * y^2)
    _test_motzkin(x, y, T, solver, set, false, true, false)
    set = PolyJuMP.SAGE.Signomials(x^2 * y)
    _test_motzkin(x, y, T, solver, set, false, false, false)
    set = PolyJuMP.SAGE.Signomials()
    _test_motzkin(x, y, T, solver, set, true, true, false)
    _test_motzkin(x, y, T, solver, set, true, false, false)
    _test_motzkin(x, y, T, solver, set, false, true, true)
    _test_motzkin(x, y, T, solver, set, false, false, true)
    set = PolyJuMP.SAGE.Polynomials(x^2 * y^2)
    _test_motzkin(x, y, T, solver, set, true, true, false)
    set = PolyJuMP.SAGE.Polynomials(x * y)
    _test_motzkin(x, y, T, solver, set, false, false, false)
    set = PolyJuMP.SAGE.Polynomials(x^4 * y^2)
    _test_motzkin(x, y, T, solver, set, false, true, false)
    set = PolyJuMP.SAGE.Polynomials(x^2 * y)
    _test_motzkin(x, y, T, solver, set, false, false, false)
    set = PolyJuMP.SAGE.Polynomials()
    _test_motzkin(x, y, T, solver, set, true, true, false)
    set = PolyJuMP.SAGE.Polynomials()
    _test_motzkin(x, y, T, solver, set, false, false, false)
    return
end

# See https://github.com/jump-dev/PolyJuMP.jl/issues/102#issuecomment-1888004697
function test_domain(x, y, T, solver)
    p = x^3 - x^2 + 2x * y - y^2 + y^3
    S = @set x >= 0 && y >= 0 && x + y >= 1
    model = Model(solver)
    setpolymodule!(model, PolyJuMP.SAGE)
    @variable(model, α)
    @objective(model, Max, α)
    @test_throws ErrorException @constraint(model, c3, p >= α, domain = S)
end

function test_optimizer_attributes(x, y, T, solver)
    # We don't specify `T` to test the fallback
    @test PolyJuMP.SAGE.Optimizer(solver) isa PolyJuMP.SAGE.Optimizer{Float64}
    optimizer = PolyJuMP.SAGE.Optimizer{T}(solver)
    @test MOI.get(optimizer, MOI.SolverName()) == "PolyJuMP.SAGE"
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
    @test MOI.get(optimizer, MOI.ResultCount()) == 0
    list = MOI.get(optimizer, MOI.Bridges.ListOfNonstandardBridges{T}())
    @test PolyJuMP.Bridges.Constraint.ToPolynomialBridge{T} in list
    @test PolyJuMP.Bridges.Objective.ToPolynomialBridge{T} in list
    @test MOI.supports_incremental_interface(optimizer)
    src = MOI.Utilities.Model{T}()
    v = MOI.add_variable(src)
    index_map = MOI.copy_to(optimizer, src)
    @test MOI.is_valid(optimizer, index_map[v])
end

function test_optimizer_motzkin(x, y, T, solver)
    model = Model(() -> PolyJuMP.SAGE.Optimizer{T}(solver))
    @variable(model, a)
    @variable(model, b)
    @objective(model, Min, a^4 * b^2 + a^2 * b^4 + 1 - 3 * a^2 * b^2)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test primal_status(model) == MOI.NO_SOLUTION
    @test result_count(model) == 0
    @test objective_bound(model) ≈ 0 atol = 1e-3
end

function test_optimizer_constrained(x, y, T, solver)
    model = Model(() -> PolyJuMP.SAGE.Optimizer{T}(solver))
    @variable(model, a)
    @objective(model, Min, a^2)
    @constraint(model, con, a^2 >= 1)
    # Setting the attribute before `optimize!` covers its `MOI.supports`
    # which is checked when the cache is copied to the optimizer
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), con, 2)
    @test MOI.get(model, PolyJuMP.MultiplierMaxdegree(), con) == 2
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ 1 rtol = 1e-3
    # Setting it after `optimize!` covers the direct forwarding to the
    # attached optimizer
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), con, 0)
    @test MOI.get(model, PolyJuMP.MultiplierMaxdegree(), con) == 0
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ 1 rtol = 1e-3
end

function test_optimizer_equality_max(x, y, T, solver)
    model = Model(() -> PolyJuMP.SAGE.Optimizer{T}(solver))
    @variable(model, a)
    @objective(model, Max, 2 - a^2)
    @constraint(model, a^2 == 1)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ 1 rtol = 1e-3
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
    TestSAGE.runtests(x, y, Float64)
end

import TypedPolynomials
@testset "DynamicPolynomials" begin
    TypedPolynomials.@polyvar(x, y)
    TestSAGE.runtests(x, y, Float64)
end
