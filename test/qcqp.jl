module TestQCQP

using Test

import MathOptInterface as MOI
import MultivariatePolynomials as MP
import PolyJuMP
import JuMP

function test_solver_name(_, _, _)
    model = Model{Float64}()
    inner = MOI.Utilities.MockOptimizer(model)
    # We don't specify `T` to test the fallback
    optimizer = PolyJuMP.QCQP.Optimizer(inner)
    @test optimizer isa PolyJuMP.QCQP.Optimizer{Float64}
    @test MOI.get(optimizer, MOI.SolverName()) == "PolyJuMP.QCQP with Mock"
end

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

MOI.Utilities.@model(
    AffineObjectiveModel,
    (),
    (MOI.LessThan, MOI.GreaterThan, MOI.EqualTo, MOI.Interval),
    (),
    (),
    (),
    (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
    (),
    (),
)

function MOI.supports(::AffineObjectiveModel{T}, ::MOI.ObjectiveFunction{F}) where {T,F<:MOI.AbstractFunction}
    return F == MOI.ScalarAffineFunction{T}
end

function MOI.supports(
    ::Model,
    ::MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction},
)
    return false
end

function _test_objective_or_constraint(x, y, T, obj::Bool)
    inner = Model{T}()
    optimizer = MOI.Utilities.MockOptimizer(inner)
    model = JuMP.GenericModel{T}(() -> PolyJuMP.QCQP.Optimizer{T}(optimizer))
    PolyJuMP.@variable(model, 1 <= a <= 2)
    PolyJuMP.@variable(model, -5 <= b <= 3)
    PolyJuMP.@constraint(model, a + b >= 1)
    if obj
        PolyJuMP.@objective(model, Min, a^3 - a^2 + 2a * b - b^2 + b^3)
    else
        PolyJuMP.@constraint(model, 0 <= a^3 - a^2 + 2a * b - b^2 + b^3 <= 1)
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
    model = JuMP.GenericModel{T}(() -> PolyJuMP.QCQP.Optimizer{T}(optimizer))
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
    @test MOI.get(inner, MOI.ConstraintFunction(), cis[3]) ≈
          o * a3 * b3 + o * a6
end

function test_no_monomials(x, y, T)
    inner = Model{T}()
    model = JuMP.GenericModel{T}() do
        return PolyJuMP.QCQP.Optimizer{T}(MOI.Utilities.MockOptimizer(inner))
    end
    PolyJuMP.@variable(model, 0 <= x[1:2] <= 2)
    PolyJuMP.@constraint(model, x[1] * x[2] == 1)
    PolyJuMP.@objective(model, Min, sum(x))
    PolyJuMP.optimize!(model)
    @test MOI.get(inner, MOI.NumberOfVariables()) == 2
    return
end

function test_scalar_constant_not_zero(x, y, T)
    inner = Model{T}()
    model = JuMP.GenericModel{T}() do
        return PolyJuMP.QCQP.Optimizer{T}(MOI.Utilities.MockOptimizer(inner))
    end
    PolyJuMP.@variable(model, -1 <= x[1:3] <= 2)
    PolyJuMP.@constraint(model, 4 * x[1] * x[2] * x[3] == 1)
    PolyJuMP.@objective(model, Min, sum(x))
    PolyJuMP.optimize!(model)
    for (F, S) in MOI.get(inner, MOI.ListOfConstraintTypesPresent())
        for ci in MOI.get(inner, MOI.ListOfConstraintIndices{F,S}())
            if F isa MOI.ScalarQuadraticFunction
                f = MOI.get(inner, MOI.ConstraintFunction(), ci)
                @test iszero(f.constant)
            end
        end
    end
    return
end

function test_unbound_polynomial(x, y, T)
    inner = Model{T}()
    model = JuMP.GenericModel{T}() do
        return PolyJuMP.QCQP.Optimizer{T}(MOI.Utilities.MockOptimizer(inner))
    end
    PolyJuMP.@variable(model, x >= 0)
    PolyJuMP.@objective(model, Min, x^3)
    PolyJuMP.optimize!(model)
    F = MOI.VariableIndex
    S = MOI.Interval{T}
    @test MOI.get(inner, MOI.NumberOfConstraints{F,S}()) == 0
    S = MOI.LessThan{T}
    @test MOI.get(inner, MOI.NumberOfConstraints{F,S}()) == 0
    S = MOI.GreaterThan{T}
    @test MOI.get(inner, MOI.NumberOfConstraints{F,S}()) == 2
    for ci in MOI.get(inner, MOI.ListOfConstraintIndices{F,S}())
        set = MOI.get(inner, MOI.ConstraintSet(), ci)
        @test set.lower == zero(T)
    end
    return
end

function test_scalar_nonlinear_function(x, y, T)
    inner = Model{T}()
    model = JuMP.GenericModel{T}() do
        return PolyJuMP.QCQP.Optimizer{T}(MOI.Utilities.MockOptimizer(inner))
    end
    PolyJuMP.@variable(model, 0 <= x <= 1)
    PolyJuMP.@expression(model, f, 0 + x)
    PolyJuMP.@expression(model, g, x^2)
    PolyJuMP.@constraint(model, f * g == 0)
    PolyJuMP.optimize!(model)
    F, S = MOI.ScalarQuadraticFunction{T}, MOI.EqualTo{T}
    @test MOI.get(inner, MOI.NumberOfConstraints{F,S}()) == 2
    @test MOI.get(inner, MOI.NumberOfVariables()) == 2
    return
end

function test_scalar_nonlinear_function_div_rem_zero(x, y, T)
    inner = Model{T}()
    model = JuMP.GenericModel{T}() do
        return PolyJuMP.QCQP.Optimizer{T}(MOI.Utilities.MockOptimizer(inner))
    end
    PolyJuMP.@variable(model, x)
    PolyJuMP.@objective(model, Min, x^3 / x)
    PolyJuMP.optimize!(model)
    @test MOI.get(inner, MOI.NumberOfVariables()) == 1
    @test isempty(MOI.get(inner, MOI.ListOfConstraintTypesPresent()))
    @test PolyJuMP.objective_function(model) isa PolyJuMP.GenericNonlinearExpr
    F = MOI.get(inner, MOI.ObjectiveFunctionType())
    @test F <: MOI.ScalarQuadraticFunction{T}
    return
end

function test_scalar_nonlinear_function_div_rem_err(x, y, T)
    inner = Model{T}()
    model = JuMP.GenericModel{T}() do
        return PolyJuMP.QCQP.Optimizer{T}(MOI.Utilities.MockOptimizer(inner))
    end
    PolyJuMP.@variable(model, x)
    PolyJuMP.@variable(model, y)
    PolyJuMP.@objective(model, Min, x^3 / y)
    @test_throws PolyJuMP.InvalidNLExpression PolyJuMP.optimize!(model)
    PolyJuMP.@objective(model, Min, 2 / y)
    @test_throws PolyJuMP.InvalidNLExpression PolyJuMP.optimize!(model)
    PolyJuMP.@objective(model, Min, 2 / (x * y^2 - y^2 * x - 1))
    PolyJuMP.optimize!(model)
    return
end

function test_scalar_nonlinear_function_div_rem_number(x, y, T)
    inner = Model{T}()
    model = JuMP.GenericModel{T}() do
        return PolyJuMP.QCQP.Optimizer{T}(MOI.Utilities.MockOptimizer(inner))
    end
    PolyJuMP.@variable(model, x)
    PolyJuMP.@objective(model, Min, 4x^3 / 2)
    PolyJuMP.optimize!(model)
    @test MOI.get(inner, MOI.NumberOfVariables()) == 2
    F, S = MOI.ScalarQuadraticFunction{T}, MOI.EqualTo{T}
    @test (F, S) in MOI.get(inner, MOI.ListOfConstraintTypesPresent())
    return
end

function test_variable_primal(x, y, T)
    inner = Model{T}()
    optimizer = MOI.Utilities.MockOptimizer(inner)
    model = JuMP.direct_generic_model(
        T,
        MOI.instantiate(
            () -> PolyJuMP.QCQP.Optimizer{T}(optimizer),
            with_bridge_type = T,
        ),
    )
    JuMP.@variable(model, 1 <= a <= 3)
    aff = JuMP.@constraint(model, a <= 1)
    cub = JuMP.@constraint(model, a^3 <= 1)
    MOI.set(model, MOI.VariablePrimal(), a, T(2))
    MOI.set(model, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(model, MOI.ConstraintDual(), aff, T(3))
    MOI.Utilities.final_touch(JuMP.backend(model), nothing)
    inner = JuMP.backend(model).model
    F = MOI.ScalarQuadraticFunction{T}
    S = MOI.LessThan{T}
    ci = first(MOI.get(inner, MOI.ListOfConstraintIndices{F,S}()))
    MOI.set(inner, MOI.ConstraintDual(), ci, T(4))
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test JuMP.value(a) == 2
    @test JuMP.dual(aff) == 3
    @test JuMP.dual(cub) == 4
end

# We test names as it's supported by `MOI.Utilities.Model`
function test_name(x, y, T)
    model = JuMP.GenericModel{T}()
    JuMP.@variable(model, 1 <= a <= 3)
    JuMP.@variable(model, 1 <= b <= 3)
    JuMP.@constraint(model, aff, a >= b)
    JuMP.@constraint(model, con_ref, a^3 >= a * b^4)
    inner = Model{T}()
    qcqp = MOI.instantiate(
        () -> PolyJuMP.QCQP.Optimizer{T}(inner),
        with_bridge_type = T,
    )
    idxmap = MOI.copy_to(qcqp, JuMP.backend(model))
    attr = MOI.VariableName()
    @test MOI.get(qcqp, attr, idxmap[JuMP.index(a)]) == "a"
    @test MOI.get(qcqp, attr, idxmap[JuMP.index(b)]) == "b"
    attr = MOI.ConstraintName()
    @test MOI.get(qcqp, attr, idxmap[JuMP.index(aff)]) == "aff"
    @test MOI.get(qcqp, attr, idxmap[JuMP.index(con_ref)]) == "con_ref"
    inner = qcqp.model.model
    F = MOI.ScalarQuadraticFunction{T}
    S = MOI.GreaterThan{T}
    ci = first(MOI.get(inner, MOI.ListOfConstraintIndices{F,S}()))
    @test_broken MOI.get(inner, attr, ci) == "con_ref"
end

function test_start(x, y, T)
    inner = MOI.Utilities.UniversalFallback(Model{T}())
    model = PolyJuMP.QCQP.Optimizer{T}(inner)
    a = MOI.add_variable(model)
    MOI.set(model, MOI.VariablePrimalStart(), a, 2one(T))
    b = MOI.add_variable(model)
    MOI.set(model, MOI.VariablePrimalStart(), b, 3one(T))
    p = PolyJuMP.ScalarPolynomialFunction(one(T) * x^3 - x * y^2, [a, b])
    ci = MOI.add_constraint(model, p, MOI.LessThan(zero(T)))
    @test MOI.is_valid(model, ci)
    MOI.Utilities.final_touch(model, nothing)
    vis = MOI.get(inner, MOI.ListOfVariableIndices())
    @test sort(MOI.get(inner, MOI.VariablePrimalStart(), vis)) == T[2, 3, 4, 9]
end

function test_inner_bridge(x, y, T)
    # The quadratic objective should be bridged after the QCQP layer
    inner = AffineObjectiveModel{T}()
    model = PolyJuMP.QCQP.Optimizer{T}(inner)
    a = MOI.add_variable(model)
    b = MOI.add_variable(model)
    p = PolyJuMP.ScalarPolynomialFunction(one(T) * x^3 - x * y^2, [a, b])
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(p)}(), p)
    MOI.Utilities.final_touch(model, nothing)
    F = MOI.ScalarAffineFunction{T}
    @test MOI.ObjectiveFunction{F}() in MOI.get(inner, MOI.ListOfModelAttributesSet())
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
