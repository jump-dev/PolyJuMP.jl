module KKT

import MutableArithmetics as MA
import MathOptInterface as MOI
import MultivariatePolynomials as MP
import SemialgebraicSets as SS
import DynamicPolynomials
import PolyJuMP

Base.@kwdef mutable struct Options{T}
    solver::Union{Nothing,SS.AbstractAlgebraicSolver} = nothing
    optimality_tolerance::Union{Nothing,T} = nothing
    feasibility_tolerance::T = Base.rtoldefault(T)
end

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    model::PolyJuMP.Model{T}
    options::Options{T}
    # Result
    solutions::Vector{PolyJuMP.Solution{T}}
    solve_time::Float64
    termination_status::MOI.TerminationStatusCode
    raw_status::String
end

function Optimizer{T}() where {T}
    optimizer = Optimizer{T}(
        PolyJuMP.Model{T}(),
        Options{T}(),
        PolyJuMP.Solution{T}[],
        NaN,
        MOI.OPTIMIZE_NOT_CALLED,
        "",
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

function MOI.supports(model::Optimizer, attr::MOI.AbstractModelAttribute)
    return MOI.supports(model.model, attr)
end

function MOI.set(model::Optimizer, attr::MOI.AbstractModelAttribute, value)
    MOI.set(model.model, attr, value)
    invalidate_solutions!(model)
    return
end

function MOI.get(
    model::Optimizer,
    attr::Union{
        MOI.AbstractModelAttribute,
        MOI.Bridges.ListOfNonstandardBridges,
    },
)
    return MOI.get(model.model, attr)
end

function MOI.supports(::Optimizer{T}, attr::MOI.RawOptimizerAttribute) where {T}
    return hasfield(Options{T}, Symbol(attr.name))
end

function MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, value)
    if !MOI.supports(model, attr)
        throw(MOI.UnsupportedAttribute(attr))
    end
    setfield!(model.options, Symbol(attr.name), value)
    return
end

function MOI.get(model::Optimizer, attr::MOI.RawOptimizerAttribute)
    if !MOI.supports(model, attr)
        throw(MOI.UnsupportedAttribute(attr))
    end
    return getfield(model.options, Symbol(attr.name))
end

function invalidate_solutions!(model::Optimizer)
    empty!(model.solutions)
    model.solve_time = NaN
    model.termination_status = MOI.OPTIMIZE_NOT_CALLED
    model.raw_status = ""
    return
end

MOI.is_valid(model::Optimizer, i::MOI.Index) = MOI.is_valid(model.model, i)

function MOI.add_variable(model::Optimizer)
    invalidate_solutions!(model)
    return MOI.add_variable(model.model)
end

function MOI.supports_constraint(
    model::Optimizer,
    ::Type{F},
    ::Type{S},
) where {F<:MOI.AbstractFunction,S<:MOI.AbstractSet}
    return MOI.supports_constraint(model.model, F, S)
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.AbstractFunction,
    set::MOI.AbstractSet,
)
    ci = MOI.add_constraint(model.model, func, set)
    invalidate_solutions!(model)
    return ci
end

MOI.supports_incremental_interface(::Optimizer) = true
function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

function _square(x::Vector{T}, n) where {T}
    return T[(i + n in eachindex(x)) ? x[i] : x[i]^2 for i in eachindex(x)]
end

function _optimize!(model::Optimizer{T}) where {T}
    lagrangian, system = PolyJuMP.lagrangian_kkt(model.model, model.options.solver)
    x = MP.variables(model.model)
    if lagrangian isa MA.Zero
        model.solutions = [
            PolyJuMP.Solution(
                zeros(T, length(x)),
                model.model,
                model.options.feasibility_tolerance,
            ),
        ]
        model.termination_status = MOI.LOCALLY_SOLVED # We could use `OPTIMAL` here but it would then make MOI tests fail as they expect `LOCALLY_SOLVED`
        model.raw_status = "Lagrangian function is zero so any solution is optimal even if the solver reports a unique solution `0`."
        return
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
            _square(sol, PolyJuMP._nineq(model.model.set)),
            model.model,
            model.options.feasibility_tolerance,
        ) for sol in solutions
    ]
    PolyJuMP.postprocess!(
        model.solutions,
        model.model,
        model.options.optimality_tolerance,
    )
    if isempty(model.solutions)
        model.termination_status = MOI.INFEASIBLE_OR_UNBOUNDED
        model.raw_status = "The KKT system is infeasible so the polynomial optimization problem is either infeasible or unbounded."
    else
        model.termination_status = MOI.LOCALLY_SOLVED
        model.raw_status = "Unless the problem is unbounded, the first of the `$(length(model.solutions))` solutions reported by the solver is optimal. The termination status is `LOCALLY_SOLVED` instead of `OPTIMAL` to leave the possibility that the problem may be unbounded."
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

MOI.get(model::Optimizer, ::MOI.ResultCount) = length(model.solutions)

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return model.solutions[attr.result_index].objective_value
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
    return model.solutions[attr.result_index].values[vi.value]
end

function _index(
    model,
    ci::MOI.ConstraintIndex{<:PolyJuMP.ScalarPolynomialFunction,<:MOI.EqualTo},
)
    return length(model.model.variables) + ci.value
end

function _index(
    model,
    ci::MOI.ConstraintIndex{
        <:PolyJuMP.ScalarPolynomialFunction,
        <:Union{MOI.LessThan,MOI.GreaterThan},
    },
)
    return length(model.model.variables) +
           SS.nequalities(model.model.set) +
           ci.value
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
    value = model.solutions[attr.result_index].values[_index(model, ci)]
    if S <: MOI.LessThan
        value = -value
    end
    return value
end

end
