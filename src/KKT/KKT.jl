module KKT

import MutableArithmetics as MA
import MathOptInterface as MOI
import MultivariatePolynomials as MP
import SemialgebraicSets as SS
import DynamicPolynomials
import PolyJuMP

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    # Model
    model::PolyJuMP.Model{T}
    # Result
    solutions::Vector{Vector{T}}
    solve_time::Float64
    termination_status::MOI.TerminationStatusCode
    raw_status::String
    # Optimizer attributes
    algebraic_solver::Union{Nothing,SS.AbstractAlgebraicSolver}
    feasibility_tolerance::T
end
function Optimizer{T}() where {T}
    optimizer = Optimizer{T}(
        PolyJuMP.Model{T}(),
        Solution{T}[],
        NaN,
        MOI.OPTIMIZE_NOT_CALLED,
        "",
        nothing,
        Base.rtoldefault(T),
        Base.rtoldefault(T),
    )
    # This can be replaced by `ListOfNonstandardBridges` once https://github.com/jump-dev/MathOptInterface.jl/issues/846 is done
    return PolyJuMP.NLToPolynomial{T}(optimizer)
end
Optimizer() = Optimizer{Float64}()

MOI.get(::Optimizer, ::MOI.SolverName) = "PolyJuMP.KKT"

MOI.is_empty(model::Optimizer) = MOI.is_empty(model.model)

function MOI.empty!(model::Optimizer)
    MOI.empty!(model.model)
    invalidate_solutions!(model)
    return
end

function MOI.support(model::Optimizer, attr::MOI.AbstractModelAttribute)
    return MOI.support(model.model, attr)
end

function MOI.set(model::Optimizer, attr::MOI.AbstractModelAttribute, value)
    MOI.set(model.model, attr, value)
    invalidate_solutions!(model)
    return
end

function MOI.get(model::Optimizer, attr::MOI.AbstractModelAttribute)
    return MOI.get(model.model, attr)
end

function MOI.set(
    model::Optimizer,
    attr::MOI.RawOptimizerAttribute,
    solver::SS.AbstractAlgebraicSolver,
)
    if attr.name != "algebraic_solver"
        throw(MOI.UnsupportedAttribute(attr))
    end
    return model.algebraic_solver = solver
end

function MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, tol)
    if attr.name == "feasibility_tolerance"
        model.feasibility_tolerance = tol
    else
        throw(MOI.UnsupportedAttribute(attr))
    end
end

function invalidate_solutions!(model::Optimizer)
    empty!(model.solution)
    model.solve_time = NaN
    model.termination_status = MOI.OPTIMIZE_NOT_CALLED
    model.raw_status = ""
    return
end

MOI.is_valid(model::Optimizer, i) = MOI.is_valid(model.model, i)

function MOI.add_variable(model::Optimizer)
    invalidate_solutions!(model)
    return MOI.add_variable(model.model)
end

function MOI.supports_constraint(model::Optimizer, ::Type{F}, ::Type{S}) where {F,S}
    return MOI.supports_constraint(model.model, F, S)
end

function MOI.add_constraint(model::Optimizer, func, set)
    ci = MOI.add_constraint(model.model, func, set)
    invalidate_solutions!(model)
    return ci
end

MOI.supports_incremental_interface(::Optimizer) = true
function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

function _add_to_system(
    system,
    lagrangian,
    ::SS.FullSpace,
    ::Bool,
)
    return lagrangian
end

function _add_to_system(
    system,
    lagrangian,
    set::SS.AlgebraicSet,
    maximization::Bool,
)
    n = SemialgebraicSets.nequalities(set)
    if iszero(n)
        return
    end
    DynamicPolynomials.@polyvar λ[1:n]
    for i in eachindex(λ)
        p = SS.equalities(set)[i]
        SS.add_equality!(system, p)
        if maximization
            lagrangian = MA.add_mul!!(lagrangian, λ[i], p)
        else
            lagrangian = MA.sub_mul!!(lagrangian, λ[i], p)
        end
    end
    return lagrangian
end

function _add_to_system(
    system,
    lagrangian,
    set::SS.BasicSemialgebraicSet,
    maximization::Bool,
)
    lagrangian = _add_to_system(system, lagrangian, set.V, maximization)
    DynamicPolynomials.@polyvar σ[1:_nineq(set)]
    for i in eachindex(σ)
        p = SS.inequalities(set)[i]
        SS.add_equality!(system, σ[i] * p)
        if maximization
            lagrangian = MA.add_mul!!(lagrangian, σ[i]^2, p)
        else
            lagrangian = MA.sub_mul!!(lagrangian, σ[i]^2, p)
        end
    end
    return lagrangian
end

function _square(x::Vector{T}, n) where {T}
    return T[(i + n in eachindex(x)) ? x[i] : x[i]^2 for i in eachindex(x)]
end

function _optimize!(model::Optimizer{T}) where {T}
    if isnothing(model.algebraic_solver)
        system = SS.AlgebraicSet{T,PolyType{T}}()
    else
        I = SS.PolynomialIdeal{T,PolyType{T}}()
        system = SS.AlgebraicSet(I, model.algebraic_solver)
    end
    if model.objective_sense == MOI.FEASIBILITY_SENSE
        lagrangian = MA.Zero()
    else
        lagrangian = MA.mutable_copy(model.objective_function)
    end
    lagrangian = _add_to_system(
        system,
        lagrangian,
        model.set,
        model.objective_sense == MOI.MAX_SENSE,
    )
    x = MP.variables(model.model)
    if lagrangian isa MA.Zero
        model.extrema = [zeros(T, length(x))]
        model.objective_values = zeros(T, 1)
        model.termination_status = MOI.LOCALLY_SOLVED # We could use `OPTIMAL` here but it would then make MOI tests fail as they expect `LOCALLY_SOLVED`
        model.raw_status = "Lagrangian function is zero so any solution is optimal even if the solver reports a unique solution `0`."
        return
    end
    ∇x = MP.differentiate(lagrangian, x)
    for p in ∇x
        SS.add_equality!(system, p)
    end
    solutions = nothing
    try # We could check `SS.is_zero_dimensional(system)` but that would only work for Gröbner basis based
        solutions = collect(system)
    catch err
        model.extrema = Vector{T}[]
        model.objective_values = zeros(T, 0)
        model.termination_status = MOI.OTHER_ERROR
        model.raw_status = "KKT system solver failed with : $(sprint(showerror, err))."
        return
    end
    model.solutions = [
        PolyJuMP.Solution(
            _square(sol, _nineq(model.set)),
            model.model,
            model.feasibility_tolerance,
        )
        for sol in solutions
    ]
    sort_unique!(solutions, model)
    if isempty(model.extrema)
        model.termination_status = MOI.INFEASIBLE_OR_UNBOUNDED
        model.raw_status = "The KKT system is infeasible so the polynomial optimization problem is either infeasible or unbounded."
    else
        model.termination_status = MOI.LOCALLY_SOLVED
        model.raw_status = "Unless the problem is unbounded, the $(length(model.extrema)) solutions reported by the solver are optimal. The termination status is `LOCALLY_SOLVED` instead of `OPTIMAL` to leave the possibility that the problem may be unbounded."
    end
    return
end

function MOI.optimize!(model::Optimizer)
    return model.solve_time = @elapsed _optimize!(model)
end

MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = model.solve_time

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    if model.termination_status === MOI.OPTIMIZE_NOT_CALLED
        return "`optimize!` has not yet been called"
    else
        return model.raw_status
    end
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    return model.termination_status
end

MOI.get(model::Optimizer, ::MOI.ResultCount) = length(model.extrema)

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return model.objective_values[attr.result_index]
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index in 1:MOI.get(model, MOI.ResultCount())
        return MOI.FEASIBLE_POINT
    else
        return MOI.NO_SOLUTION
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.throw_if_not_valid(model, vi)
    MOI.check_result_index_bounds(model, attr)
    return model.extrema[attr.result_index][vi.value]
end

function _index(
    model,
    ci::MOI.ConstraintIndex{<:PolyJuMP.ScalarPolynomialFunction,<:MOI.EqualTo},
)
    return length(model.variables) + ci.value
end

function _index(
    model,
    ci::MOI.ConstraintIndex{
        <:PolyJuMP.ScalarPolynomialFunction,
        <:Union{MOI.LessThan,MOI.GreaterThan},
    },
)
    return length(model.variables) + SS.nequalities(model.set) + ci.value
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if attr.result_index in 1:MOI.get(model, MOI.ResultCount())
        return MOI.FEASIBLE_POINT
    else
        return MOI.NO_SOLUTION
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{F,S},
) where {F,S}
    MOI.throw_if_not_valid(model, ci)
    MOI.check_result_index_bounds(model, attr)
    value = model.extrema[attr.result_index][_index(model, ci)]
    if S <: MOI.LessThan
        value = -value
    end
    return value
end

end
