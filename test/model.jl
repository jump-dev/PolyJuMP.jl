module TestModel

using Test

import LinearAlgebra
import MultivariatePolynomials as MP
import ECOS
using JuMP
using PolyJuMP

# Polynomial matrix `[1 x^2; x^2 1]` as `MOI.VectorQuadraticFunction` with
# the variable `x` of index `1`
function _quadratic_matrix()
    model = JuMP.Model()
    JuMP.@variable(model, x)
    return JuMP.moi_function([1, x^2, 1])
end

function test_matrix_inequality_triangle()
    model = PolyJuMP.Model{Float64}()
    @test MOI.supports_constraint(
        model,
        MOI.VectorQuadraticFunction{Float64},
        MOI.PositiveSemidefiniteConeTriangle,
    )
    MOI.add_variable(model)
    set = MOI.PositiveSemidefiniteConeTriangle(2)
    ci = MOI.add_constraint(model, _quadratic_matrix(), set)
    @test MOI.is_valid(model, ci)
    @test !MOI.is_empty(model)
    @test length(model.matrix_inequalities) == 1
    x = only(MP.variables(model))
    @test model.matrix_inequalities[1] == [1 x^2; x^2 1]
    @test model.set isa PolyJuMP.SS.FullSpace
    # The minimum eigenvalue of `[1 x^2; x^2 1]` is `1 - x^2`
    sol = PolyJuMP.Solution([0.5], model, 1e-6)
    @test iszero(sol.max_constraint_violation)
    @test sol.status == MOI.FEASIBLE_POINT
    sol = PolyJuMP.Solution([2.0], model, 1e-6)
    @test sol.max_constraint_violation ≈ 3
    @test sol.status == MOI.INFEASIBLE_POINT
    MOI.empty!(model)
    @test MOI.is_empty(model)
    @test !MOI.is_valid(model, ci)
    return
end

function test_matrix_inequality_nonlinear()
    model = PolyJuMP.Model{Float64}()
    vi = MOI.add_variable(model)
    func = MOI.VectorNonlinearFunction([
        MOI.ScalarNonlinearFunction(:+, Any[1.0]),
        MOI.ScalarNonlinearFunction(:^, Any[vi, 2]),
        MOI.ScalarNonlinearFunction(:+, Any[1.0]),
    ])
    set = MOI.PositiveSemidefiniteConeTriangle(2)
    @test MOI.supports_constraint(model, typeof(func), typeof(set))
    ci = MOI.add_constraint(model, func, set)
    @test MOI.is_valid(model, ci)
    x = only(MP.variables(model))
    @test only(model.matrix_inequalities) == [1 x^2; x^2 1]
    return
end

include("testpolymodule.jl")

# Relaxation optimizer supporting matrix inequalities with the cones of
# `DummyPolyModule`. As these cones are not supported by any solver, the
# relaxation can only be solved by a `MOI.Utilities.MockOptimizer`.
mutable struct MatrixOptimizer{T} <: PolyJuMP.AbstractRelaxationOptimizer{T}
    model::PolyJuMP.Model{T}
    multiplier_maxdegree::Dict{MOI.ConstraintIndex,Int}
    solver::Any
    relaxation::Union{Nothing,JuMP.GenericModel{T}}
    solutions::Vector{PolyJuMP.Solution{T}}
    feasibility_tolerance::T
    solve_time::Float64
end

function MatrixOptimizer{T}(solver = nothing) where {T}
    return MatrixOptimizer{T}(
        PolyJuMP.Model{T}(),
        Dict{MOI.ConstraintIndex,Int}(),
        solver,
        nothing,
        PolyJuMP.Solution{T}[],
        sqrt(eps(T)),
        NaN,
    )
end

PolyJuMP.nonnegativity_cone(::MatrixOptimizer) = DummyPolyModule.DummyNonNeg()

function PolyJuMP.matrix_nonnegativity_cone(::MatrixOptimizer)
    return DummyPolyModule.DummyPosDefMatrix()
end

function PolyJuMP.recover_solutions(
    ::MatrixOptimizer{T},
    relaxation,
    cref,
    lagrangian,
) where {T}
    return PolyJuMP.Solution{T}[]
end

function _mock_solver()
    return MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
end

function test_matrix_inequality_relaxation()
    model = JuMP.Model(() -> MatrixOptimizer{Float64}(_mock_solver))
    JuMP.@variable(model, x)
    matrix = JuMP.@constraint(
        model,
        LinearAlgebra.Symmetric([1 x; x 1]) in JuMP.PSDCone(),
    )
    JuMP.@objective(model, Min, x)
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), matrix, 2)
    @test MOI.get(model, PolyJuMP.MultiplierMaxdegree(), matrix) == 2
    JuMP.optimize!(model)
    relaxation = JuMP.unsafe_backend(model).relaxation
    F, S = only(filter(JuMP.list_of_constraint_types(relaxation)) do (_, S)
        return S <: DummyPolyModule.PosDefMatrix
    end)
    cons = JuMP.all_constraints(relaxation, F, S)
    # The multiplier of the matrix inequality is a `2 × 2` matrix of
    # polynomials of degree `2`
    set = JuMP.constraint_object(only(cons)).set
    @test size(set.monomials) == (2, 2)
    @test all(m -> MP.maxdegree(m) == 2, set.monomials)
    return
end

function test_multiplier_maxdegree_square_bridge()
    optimizer = MatrixOptimizer{Float64}()
    model = MOI.Bridges.Constraint.Square{Float64}(optimizer)
    vi = MOI.add_variable(model)
    F = MOI.VectorAffineFunction{Float64}
    S = MOI.PositiveSemidefiniteConeSquare
    @test MOI.supports_constraint(model, F, S)
    @test MOI.supports(
        model,
        PolyJuMP.MultiplierMaxdegree(),
        MOI.ConstraintIndex{F,S},
    )
    func = MOI.Utilities.operate(vcat, Float64, 1.0, 1.0 * vi, 1.0 * vi, 1.0)
    ci = MOI.add_constraint(model, func, S(2))
    @test MOI.Bridges.is_bridged(model, ci)
    @test isnothing(MOI.get(model, PolyJuMP.MultiplierMaxdegree(), ci))
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), ci, 2)
    @test MOI.get(model, PolyJuMP.MultiplierMaxdegree(), ci) == 2
    @test only(values(optimizer.multiplier_maxdegree)) == 2
    x = only(MP.variables(optimizer.model))
    @test only(optimizer.model.matrix_inequalities) == [1 x; x 1]
    return
end

function test_matrix_inequality_not_supported()
    optimizer = PolyJuMP.SAGE.Optimizer(ECOS.Optimizer)
    @test isnothing(PolyJuMP.matrix_nonnegativity_cone(optimizer))
    F = MOI.VectorQuadraticFunction{Float64}
    S = MOI.PositiveSemidefiniteConeTriangle
    @test !MOI.supports_constraint(optimizer, F, S)
    @test MOI.get(optimizer, MOI.ConstraintBridgingCost{F,S}()) == Inf
    return
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestModel.runtests()
