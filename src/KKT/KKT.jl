module KKT

import MutableArithmetics
const MA = MutableArithmetics

import MathOptInterface
const MOI = MathOptInterface

import MultivariatePolynomials
const MP = MultivariatePolynomials

import SemialgebraicSets

import DynamicPolynomials

import PolyJuMP

const VariableOrder =
    DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}
const MonomialOrder = DynamicPolynomials.Graded{MP.LexOrder}
const VarType = DynamicPolynomials.Variable{VariableOrder,MonomialOrder}
const PolyType{T} = DynamicPolynomials.Polynomial{VariableOrder,MonomialOrder,T}

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    # Model
    variables::Dict{MOI.VariableIndex,VarType}
    objective_sense::MOI.OptimizationSense
    objective_function::Union{Nothing,PolyType{T}}
    set::Any
    # Result
    extrema::Vector{Vector{T}}
    objective_values::Vector{T}
    solve_time::Float64
    termination_status::MOI.TerminationStatusCode
    raw_status::String
    # Optimizer attributes
    algebraic_solver::Union{Nothing,SemialgebraicSets.AbstractAlgebraicSolver}
    feasibility_tolerance::T
    optimality_tolerance::T
end
function Optimizer{T}() where {T}
    optimizer = Optimizer{T}(
        Dict{MOI.VariableIndex,VarType}(),
        MOI.FEASIBILITY_SENSE,
        nothing,
        SemialgebraicSets.FullSpace(),
        Vector{T}[],
        T[],
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

function MOI.get(
    ::Optimizer{T},
    ::MOI.Bridges.ListOfNonstandardBridges{T},
) where {T}
    return [
        PolyJuMP.Bridges.Constraint.ToPolynomialBridge{T},
        PolyJuMP.Bridges.Objective.ToPolynomialBridge{T},
    ]
end

function MOI.set(
    model::Optimizer,
    attr::MOI.RawOptimizerAttribute,
    solver::SemialgebraicSets.AbstractAlgebraicSolver,
)
    if attr.name != "algebraic_solver"
        throw(MOI.UnsupportedAttribute(attr))
    end
    return model.algebraic_solver = solver
end

function MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, tol)
    if attr.name == "feasibility_tolerance"
        model.feasibility_tolerance = tol
    elseif attr.name != "optimality_tolerance"
        model.optimality_tolerance = tol
    else
        throw(MOI.UnsupportedAttribute(attr))
    end
end

function MOI.is_empty(model::Optimizer)
    return isempty(model.variables) &&
           model.objective_sense == MOI.FEASIBILITY_SENSE &&
           isnothing(model.objective_function) &&
           model.set isa SemialgebraicSets.FullSpace
end

function MOI.empty!(model::Optimizer)
    empty!(model.variables)
    model.objective_sense = MOI.FEASIBILITY_SENSE
    model.objective_function = nothing
    model.set = SemialgebraicSets.FullSpace()
    empty!(model.extrema)
    empty!(model.objective_values)
    model.solve_time = NaN
    model.termination_status = MOI.OPTIMIZE_NOT_CALLED
    model.raw_status = ""
    return
end

function MOI.is_valid(model::Optimizer, vi::MOI.VariableIndex)
    return in(vi.value, 1:length(model.variables))
end

function MOI.add_variable(model::Optimizer)
    i = length(model.variables) + 1
    vi = MOI.VariableIndex(i)
    var = DynamicPolynomials.Variable("x[$i]", VariableOrder, MonomialOrder)
    model.variables[vi] = var
    model.termination_status = MOI.OPTIMIZE_NOT_CALLED
    return vi
end

function _polynomial(variables, func::PolyJuMP.ScalarPolynomialFunction)
    new_variables = VarType[variables[vi] for vi in func.variables]
    return func.polynomial(MP.variables(func.polynomial) => new_variables)
end

function _set(poly::MP.AbstractPolynomialLike, set::MOI.EqualTo)
    return SemialgebraicSets.equality(poly, MOI.constant(set))
end

function _set(poly::MP.AbstractPolynomialLike, set::MOI.GreaterThan)
    return SemialgebraicSets.PolynomialInequality(poly - MOI.constant(set))
end

function _set(poly::MP.AbstractPolynomialLike, set::MOI.LessThan)
    return SemialgebraicSets.PolynomialInequality(MOI.constant(set) - poly)
end

_nineq(::SemialgebraicSets.AbstractAlgebraicSet) = 0
_nineq(set) = SemialgebraicSets.ninequalities(set)

_num(set, ::Type{<:MOI.EqualTo}) = SemialgebraicSets.nequalities(set)
function _num(set, ::Type{<:Union{MOI.LessThan,MOI.GreaterThan}})
    return _nineq(set)
end

function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{<:PolyJuMP.ScalarPolynomialFunction{T}},
    ::Type{<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T}}},
) where {T}
    return true
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{<:PolyJuMP.ScalarPolynomialFunction,S},
) where {S}
    return ci.value in 1:_num(model.set, S)
end

function MOI.add_constraint(
    model::Optimizer{T},
    func::PolyJuMP.ScalarPolynomialFunction{T},
    set::MOI.AbstractScalarSet,
) where {T}
    model.set = model.set ∩ _set(_polynomial(model.variables, func), set)
    i = _num(model.set, typeof(set))
    model.termination_status = MOI.OPTIMIZE_NOT_CALLED
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(i)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    model.objective_sense = sense
    model.termination_status = MOI.OPTIMIZE_NOT_CALLED
    return
end

function MOI.supports(
    ::Optimizer{T},
    ::MOI.ObjectiveFunction{<:PolyJuMP.ScalarPolynomialFunction{T}},
) where {T}
    return true
end
function MOI.set(
    model::Optimizer{T},
    ::MOI.ObjectiveFunction,
    func::PolyJuMP.ScalarPolynomialFunction{T},
) where {T}
    model.objective_function = _polynomial(model.variables, func)
    model.termination_status = MOI.OPTIMIZE_NOT_CALLED
    return
end

function _add_to_system(
    system,
    lagrangian,
    set::SemialgebraicSets.FullSpace,
    ::Bool,
)
    return lagrangian
end

function _add_to_system(
    system,
    lagrangian,
    set::SemialgebraicSets.AlgebraicSet,
    maximization::Bool,
)
    n = SemialgebraicSets.nequalities(set)
    if iszero(n)
        return
    end
    DynamicPolynomials.@polyvar λ[1:n]
    for i in eachindex(λ)
        p = SemialgebraicSets.equalities(set)[i]
        SemialgebraicSets.add_equality!(system, p)
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
    set::SemialgebraicSets.BasicSemialgebraicSet,
    maximization::Bool,
)
    lagrangian = _add_to_system(system, lagrangian, set.V, maximization)
    DynamicPolynomials.@polyvar σ[1:_nineq(set)]
    for i in eachindex(σ)
        p = SemialgebraicSets.inequalities(set)[i]
        SemialgebraicSets.add_equality!(system, σ[i] * p)
        if maximization
            lagrangian = MA.add_mul!!(lagrangian, σ[i]^2, p)
        else
            lagrangian = MA.sub_mul!!(lagrangian, σ[i]^2, p)
        end
    end
    return lagrangian
end

function _violates_inequalities(
    set::Union{SemialgebraicSets.FullSpace,SemialgebraicSets.AlgebraicSet},
    x,
    sol,
    tol,
)
    return false
end

function _violates_inequalities(
    set::SemialgebraicSets.BasicSemialgebraicSet,
    x,
    sol,
    tol,
)
    return any(p -> p(x => sol) < -tol, SemialgebraicSets.inequalities(set))
end

function _square(x::Vector{T}, n) where {T}
    return T[(i + n in eachindex(x)) ? x[i] : x[i]^2 for i in eachindex(x)]
end

MOI.supports_incremental_interface(::Optimizer) = true
function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

function _optimize!(model::Optimizer{T}) where {T}
    if isnothing(model.algebraic_solver)
        system = SemialgebraicSets.AlgebraicSet{T,PolyType{T}}()
    else
        I = SemialgebraicSets.PolynomialIdeal{T,PolyType{T}}()
        system = SemialgebraicSets.AlgebraicSet(I, model.algebraic_solver)
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
    x = sort!(collect(values(model.variables)), rev = true)
    if lagrangian isa MA.Zero
        model.extrema = [zeros(T, length(x))]
        model.objective_values = zeros(T, 1)
        model.termination_status = MOI.LOCALLY_SOLVED # We could use `OPTIMAL` here but it would then make MOI tests fail as they expect `LOCALLY_SOLVED`
        model.raw_status = "Lagrangian function is zero so any solution is optimal even if the solver reports a unique solution `0`."
        return
    end
    ∇x = MP.differentiate(lagrangian, x)
    for p in ∇x
        SemialgebraicSets.add_equality!(system, p)
    end
    solutions = nothing
    try # We could check `SemialgebraicSets.is_zero_dimensional(system)` but that would only work for Gröbner basis based
        solutions = collect(system)
    catch err
        model.extrema = Vector{T}[]
        model.objective_values = zeros(T, 0)
        model.termination_status = MOI.OTHER_ERROR
        model.raw_status = "KKT system solver failed with : $(sprint(showerror, err))."
        return
    end
    model.extrema = Vector{T}[
        _square(sol, _nineq(model.set)) for
        sol in solutions if !_violates_inequalities(
            model.set,
            x,
            sol,
            model.feasibility_tolerance,
        )
    ]
    if model.objective_sense != MOI.FEASIBILITY_SENSE
        model.objective_values = T[
            model.objective_function(x => sol[eachindex(x)]) for
            sol in model.extrema
        ]
        I = sortperm(
            model.objective_values,
            rev = model.objective_sense == MOI.MAX_SENSE,
        )
        model.extrema = model.extrema[I]
        model.objective_values = model.objective_values[I]
        # Even if SemialgebraicSets remove duplicates, we may have solution with different `σ` but same `σ^2`
        J = Int[
            i for i in eachindex(model.extrema) if
            i == 1 || !isapprox(model.extrema[i], model.extrema[i-1])
        ]
        model.extrema = model.extrema[J]
        model.objective_values = model.objective_values[J]
        i = findfirst(eachindex(model.objective_values)) do i
            return abs(model.objective_values[i] - first(model.objective_values)) > model.optimality_tolerance
        end
        if !isnothing(i)
            model.extrema = model.extrema[1:(i-1)]
            model.objective_values = model.objective_values[1:(i-1)]
        end
    end
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
    return length(model.variables) +
           SemialgebraicSets.nequalities(model.set) +
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
    value = model.extrema[attr.result_index][_index(model, ci)]
    if S <: MOI.LessThan
        value = -value
    end
    return value
end

end
