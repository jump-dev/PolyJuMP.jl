module TestToPolynomial

using Test
import DynamicPolynomials
import MathOptInterface as MOI
import MultivariatePolynomials as MP
import PolyJuMP

const CB = PolyJuMP.Bridges.Constraint
const OB = PolyJuMP.Bridges.Objective

function _native_functions(::Type{T}, x, y) where {T}
    DynamicPolynomials.@polyvar p q
    return (
        x,
        T(2) * x + T(3) * y + T(1),
        T(2) * x * x + T(3) * y + T(1),
        PolyJuMP.ScalarPolynomialFunction(
            MP.polynomial(T(2) * p^3 + T(3) * q + T(1)),
            [x, y],
        ),
    )
end

function _nonlinear_function(::Type{T}, x, y) where {T}
    return MOI.ScalarNonlinearFunction(
        :+,
        Any[
            MOI.ScalarNonlinearFunction(:^, Any[x, 3]),
            MOI.ScalarNonlinearFunction(:*, Any[T(2), y]),
            T(1),
        ],
    )
end

function _test_unchanged(actual, expected)
    @test typeof(actual) === typeof(expected)
    if expected isa PolyJuMP.ScalarPolynomialFunction
        @test actual.variables == expected.variables
        @test actual.polynomial == expected.polynomial
    elseif expected isa MOI.VariableIndex
        @test actual == expected
    else
        @test actual ≈ expected
    end
    return
end

function _test_polynomial(f, x, y, ::Type{T}) where {T}
    @test f isa PolyJuMP.ScalarPolynomialFunction{T}
    @test Set(f.variables) == Set([x, y])
    variables = Dict(zip(f.variables, MP.variables(f.polynomial)))
    @test f.polynomial == variables[x]^3 + 2variables[y] + 1
    return
end

function test_constraints(::Type{T}) where {T}
    for set in (MOI.LessThan(T(5)), MOI.EqualTo(T(5)))
        backend = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}())
        model = MOI.Bridges.Constraint.SingleBridgeOptimizer{
            CB.ToPolynomialBridge{T},
        }(
            backend,
        )
        x, y = MOI.add_variables(model, 2)
        for f in _native_functions(T, x, y)
            @test !MOI.supports_constraint(
                CB.ToPolynomialBridge{T},
                typeof(f),
                typeof(set),
            )
            ci = MOI.add_constraint(model, f, set)
            @test !MOI.Bridges.is_bridged(model, ci)
            @test MOI.is_valid(backend, ci)
            _test_unchanged(MOI.get(backend, MOI.ConstraintFunction(), ci), f)
            @test MOI.get(backend, MOI.ConstraintSet(), ci) == set
        end
        f = _nonlinear_function(T, x, y)
        @test MOI.supports_constraint(
            CB.ToPolynomialBridge{T},
            typeof(f),
            typeof(set),
        )
        ci = MOI.add_constraint(model, f, set)
        @test MOI.Bridges.is_bridged(model, ci)
        bridge = MOI.Bridges.bridge(model, ci)
        attr = MOI.ListOfConstraintIndices{CB.FuncType{T},typeof(set)}()
        target = only(MOI.get(bridge, attr))
        _test_polynomial(
            MOI.get(backend, MOI.ConstraintFunction(), target),
            x,
            y,
            T,
        )
        @test MOI.get(backend, MOI.ConstraintSet(), target) == set
    end
    return
end

function test_objectives(::Type{T}) where {T}
    backend = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}())
    model =
        MOI.Bridges.Objective.SingleBridgeOptimizer{OB.ToPolynomialBridge{T}}(
            backend,
        )
    x, y = MOI.add_variables(model, 2)
    functions = _native_functions(T, x, y)
    for f in functions
        @test !MOI.Bridges.Objective.supports_objective_function(
            OB.ToPolynomialBridge{T},
            typeof(f),
        )
        @test !MOI.Bridges.is_bridged(model, typeof(f))
        attr = MOI.ObjectiveFunction{typeof(f)}()
        MOI.set(model, attr, f)
        @test MOI.get(backend, MOI.ObjectiveFunctionType()) === typeof(f)
        _test_unchanged(MOI.get(backend, attr), f)
    end
    for sense in (MOI.MIN_SENSE, MOI.MAX_SENSE)
        f = _nonlinear_function(T, x, y)
        @test MOI.Bridges.Objective.supports_objective_function(
            OB.ToPolynomialBridge{T},
            typeof(f),
        )
        MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
        MOI.set(model, MOI.ObjectiveSense(), sense)
        @test MOI.get(model, MOI.ObjectiveFunctionType()) === typeof(f)
        @test MOI.get(backend, MOI.ObjectiveFunctionType()) === OB.FuncType{T}
        @test MOI.get(backend, MOI.ObjectiveSense()) == sense
        _test_polynomial(
            MOI.get(backend, MOI.ObjectiveFunction{OB.FuncType{T}}()),
            x,
            y,
            T,
        )
        # Switching back to an affine objective must remove the bridge.
        affine = functions[2]
        attr = MOI.ObjectiveFunction{typeof(affine)}()
        MOI.set(model, attr, affine)
        @test MOI.get(model, MOI.ObjectiveFunctionType()) === typeof(affine)
        @test MOI.get(backend, MOI.ObjectiveFunctionType()) === typeof(affine)
        _test_unchanged(MOI.get(backend, attr), affine)
    end
    return
end

@testset "ToPolynomialBridge $T" for T in (Int, Float64)
    @testset "Constraints" begin
        test_constraints(T)
    end
    @testset "Objectives" begin
        test_objectives(T)
    end
end

end
