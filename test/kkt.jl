module TestKKT

using Test

using MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets

import DynamicPolynomials

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
    @test _dual(model, c1) ≈ -√2
    @test _dual(model, c2) ≈ -√2 / 2
end

function test_algebraic(var, T, solver)
    model = PolyJuMP.KKT.Optimizer{T}()
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), solver)
    t = MP.similar_variable(var, Val{:t})
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
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
        PolyJuMP.ScalarPolynomialFunction(o * x^2 + o * y^2 - o * t^(2), vars),
        MOI.LessThan(z),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj = PolyJuMP.ScalarPolynomialFunction(o * x + o * y, vars[2:3])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    return _test_solution(model, vars, c1, c2)
end

function _test_linquad(T, F1, O, solver)
    model = MOI.instantiate(
        optimizer_with_attributes(
            PolyJuMP.KKT.Optimizer{T},
            "solver" => solver,
        ),
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

function test_linquad(var, T, solver)
    if !(var isa DynamicPolynomials.Variable)
        return # Avoid running the same thing several times
    end
    for F1 in [MOI.VariableIndex, MOI.ScalarAffineFunction]
        for O in [MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction]
            _test_linquad(T, F1, O, solver)
        end
    end
end

function _test_JuMP(F1, O, solver)
    model = Model(PolyJuMP.KKT.Optimizer)
    set_optimizer_attribute(model, "solver", solver)
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
        @NLobjective(model, Max, 1 * x + 1 * y^(1.0))
    end
    return _test_solution(model, [t, x, y], c1, c2)
end

function test_JuMP(var, T, solver)
    if !(var isa DynamicPolynomials.Variable) || T != Float64
        return # Avoid running the same thing several times
    end
    for F1 in [MOI.VariableIndex, MOI.ScalarAffineFunction, 1, 2]
        for O in
            [MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction, 1, 2, 3, 4]
            _test_JuMP(F1, O, solver)
        end
    end
end

function test_MOI_runtests(var, T, solver)
    if !(var isa DynamicPolynomials.Variable) || T != Float64
        return # Avoid running the same thing several times
    end
    config = MOI.Test.Config(
        rtol = 1e-6,
        atol = 1e-6,
        optimal_status = MOI.LOCALLY_SOLVED,
        exclude = Any[MOI.SolverVersion, MOI.ObjectiveBound],
    )
    optimizer = MOI.instantiate(PolyJuMP.KKT.Optimizer{T}, with_bridge_type = T)
    @test MOI.get(optimizer, MOI.SolverName()) == "PolyJuMP.KKT"
    MOI.set(optimizer, MOI.RawOptimizerAttribute("solver"), solver)
    # Remove `ZerosBridge` otherwise querying `ConstraintDual` won't work in `test_quadratic_constraint_GreaterThan`
    MOI.Bridges.remove_bridge(optimizer, MOI.Bridges.Variable.ZerosBridge{T})
    cache = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}())
    cached = MOI.Utilities.CachingOptimizer(cache, optimizer)
    exclude = [
        # See https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/issues/522
        "test_solve_TerminationStatus_DUAL_INFEASIBLE",
        ### Non zero dimensional KKT system
        ## Infinite set of optimal solution:
        # min_x 0 | x ≥ 1
        "test_objective_FEASIBILITY_SENSE_clears_objective",
        # min 2x^2 | x ≥ 1, y ≥ 2
        "test_objective_qp_ObjectiveFunction_edge_cases",
        ## Feasibility sense
        # When there is no objective, we can just take all duals to be zero and then
        # the primal just need to satisfy the equality constraints.
        "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_EqualTo_lower",
        "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_EqualTo_upper",
        "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_GreaterThan",
        "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_lower",
        "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_LessThan",
        "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_upper",
        #"test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_LessThan",
        r"test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_VariableIndex_LessThan$",
        r"test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_VariableIndex_LessThan_max$",
        "test_linear_FEASIBILITY_SENSE",
        ## Misc
        # min t
        # s.t. x + y >= 1 (c1)
        #      x^2 + y^2 <= t^2 (c2)
        #      t >= 0 (bound)
        # `c1 = c2 = t = 0`, `bound = ±1`, `x, y` can be anything so the KKT system is not zero-dimensional
        "test_quadratic_SecondOrderCone_basic",
        # Max x  + y
        # st -x  + y >= 0 (c1[1])
        #     x  + y >= 0 (c1[2])
        #     x² + y <= 2 (c2)
        # `c[1]^2 = -1`, `c1[1] = 0 = c2`, then `x, y` can be anything as long as `x = -y`.
        "test_quadratic_constraint_integration",
        # min x
        # s.t. x >= 1 (xl)
        #      x <= 1 (xu)
        # `x = 1` and `xl, xu` can be anything as long as `xl - xu = 1`
        "test_solve_VariableIndex_ConstraintDual_MAX_SENSE",
        "test_solve_VariableIndex_ConstraintDual_MIN_SENSE",
        # `max x + 2y | y + x^2 + x^2 <= 1, x >= 0.5, y >= 0.5`
        # With `x = y = 0.5`, the three constraints are tight so the gradient has only 2 equations for 3 multipliers -> not zero-dimensional
        "test_constraint_qcp_duplicate_diagonal",
        # Not algebraic
        "test_nonlinear_expression_hs110",
        "test_nonlinear_expression_hs109",
        # System too big for both Buchberger and HomotopyContinuation
        r"test_nonlinear_duals$",
        r"test_nonlinear_expression_hs071$",
        r"test_nonlinear_expression_hs071_epigraph$",
    ]
    if solver isa typeof(SemialgebraicSets.default_algebraic_solver(Float64))
    else
        # FIXME
        # test_nonlinear_expression_quartic: Test Failed at /home/blegat/.julia/packages/MathOptInterface/BlCD1/src/Test/test_nonlinear.jl:1378
        #   Expression: MOI.get(model, MOI.TerminationStatus()) == config.optimal_status
        #    Evaluated: MathOptInterface.INFEASIBLE_OR_UNBOUNDED == MathOptInterface.LOCALLY_SOLVED
        push!(exclude, "test_nonlinear_expression_quartic")
    end
    MOI.Test.runtests(cached, config; exclude)
    return
end

const SOLVERS = Any[SemialgebraicSets.default_algebraic_solver(Float64),]

@static if Sys.WORD_SIZE == 64 # Issue with 32 bits, see https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/issues/476
    import HomotopyContinuation
    push!(
        SOLVERS,
        HomotopyContinuation.SemialgebraicSetsHCSolver(; compile = false),
    )
end

function runtests(var, T)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name) $(typeof(solver))" for solver in SOLVERS
                getfield(@__MODULE__, name)(var, T, solver)
            end
        end
    end
end

end # module

import DynamicPolynomials
DynamicPolynomials.@polyvar(x)
TestKKT.runtests(x, Float64)
import TypedPolynomials
TypedPolynomials.@polyvar(y)
TestKKT.runtests(y, Float64)
