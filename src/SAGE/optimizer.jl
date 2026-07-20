"""
    Optimizer{T}(solver)

Optimizer computing a bound on the objective value of a polynomial
optimization problem using its SAGE relaxation, see
[`PolyJuMP.AbstractRelaxationOptimizer`](@ref) with the cone
[`Polynomials`](@ref) as certificate of nonnegativity.
The relaxation is solved with `solver` and the bound can be queried with
`MOI.ObjectiveBound()`.
"""
mutable struct Optimizer{T} <: PolyJuMP.AbstractRelaxationOptimizer{T}
    model::PolyJuMP.Model{T}
    multiplier_maxdegree::Dict{MOI.ConstraintIndex,Int}
    solver::Any
    relaxation::Union{Nothing,JuMP.GenericModel{T}}
    solve_time::Float64
end

function Optimizer{T}(solver) where {T}
    return Optimizer{T}(
        PolyJuMP.Model{T}(),
        Dict{MOI.ConstraintIndex,Int}(),
        solver,
        nothing,
        NaN,
    )
end

Optimizer(solver) = Optimizer{Float64}(solver)

MOI.get(::Optimizer, ::MOI.SolverName) = "PolyJuMP.SAGE"

PolyJuMP.nonnegativity_cone(::Optimizer) = Polynomials()
