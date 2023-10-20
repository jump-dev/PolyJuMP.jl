mutable struct Model{T} <: MOI.ModelLike
    variables::Dict{MOI.VariableIndex,VarType}
    objective_sense::MOI.OptimizationSense
    objective_function::Union{Nothing,PolyType{T}}
    set::Any
    function Model{T}() where {T}
        return new{T}(
            Dict{MOI.VariableIndex,VarType}(),
            MOI.FEASIBILITY_SENSE,
            nothing,
            SS.FullSpace(),
        )
    end
end

function MP.variables(model::Model)
    return sort!(collect(values(model.variables)), rev = true)
end

struct Solution{T}
    values::Vector{T}
    objective_value::T
    max_constraint_violation::T
    status::MOI.ResultStatusCode
end

function Solution(
    values::Vector{T},
    model::Model{T},
    feasibility_tolerance::T,
) where {T}
    ε = _max_constraint_violation(model.set, MP.variables(model), values)
    x = MP.variables(model)
    obj = if isnothing(model.objective_function)
        zero(T)
    else
        model.objective_function(x => values[eachindex(x)])
    end
    status = if ε < feasibility_tolerance
        MOI.FEASIBLE_POINT
    elseif ε < 100feasibility_tolerance
        MOI.NEARLY_FEASIBLE_POINT
    else
        MOI.INFEASIBLE_POINT
    end
    return Solution{T}(values, obj, ε, status)
end

function _max_constraint_violation(
    ::SS.FullSpace,
    x,
    sol::AbstractVector{T},
) where {T}
    return zero(T)
end

function _max_constraint_violation(set::SS.AlgebraicSet, x, sol)
    return maximum(p -> abs(p(x => sol)), SS.equalities(set))
end

function _max_constraint_violation(set::SS.BasicSemialgebraicSet, x, sol)
    return max(
        _max_constraint_violation(set.V, x, sol),
        maximum(p -> -p(x => sol), SS.inequalities(set)),
    )
end

function _status_priority(r::MOI.ResultStatusCode)
    if r == MOI.FEASIBLE_POINT
        return 0
    elseif r == MOI.NEARLY_FEASIBLE_POINT
        return 1
    else
        @assert r == MOI.INFEASIBLE_POINT
        return 2
    end
end

function _priority(r::Solution, sense::MOI.OptimizationSense)
    obj = r.objective_value
    if sense == MOI.MAX_SENSE
        obj = -obj
    end
    return (_status_priority(r.status), obj)
end

function postprocess!(
    solutions::Vector{<:Solution},
    model::Model,
    optimality_tolerance,
)
    sort!(solutions; by = Base.Fix2(_priority, model.objective_sense))
    # Even if SemialgebraicSets remove duplicates, we may have solution with different `σ` but same `σ^2`
    J = Int[
        i for i in eachindex(solutions) if
        i != 1 && isapprox(solutions[i].values, solutions[i-1].values)
    ]
    deleteat!(solutions, J)
    if !isempty(solutions) && !isnothing(optimality_tolerance)
        sign = model.objective_sense == MOI.MAX_SENSE ? -1 : 1
        best = first(solution).objective_value
        filter!(solutions) do sol
            return sign(sol.objective_value - best) <= optimality_tolerance
        end
    end
    return
end

function MOI.get(
    ::Model{T},
    ::MOI.Bridges.ListOfNonstandardBridges{T},
) where {T}
    return [
        Bridges.Constraint.ToPolynomialBridge{T},
        Bridges.Objective.ToPolynomialBridge{T},
    ]
end

function MOI.is_empty(model::Model)
    return isempty(model.variables) &&
           model.objective_sense == MOI.FEASIBILITY_SENSE &&
           isnothing(model.objective_function) &&
           model.set isa SS.FullSpace
end

function MOI.empty!(model::Model)
    empty!(model.variables)
    model.objective_sense = MOI.FEASIBILITY_SENSE
    model.objective_function = nothing
    model.set = SS.FullSpace()
    return
end

function MOI.is_valid(model::Model, vi::MOI.VariableIndex)
    return in(vi.value, 1:length(model.variables))
end

function MOI.add_variable(model::Model)
    i = length(model.variables) + 1
    vi = MOI.VariableIndex(i)
    var = DynamicPolynomials.Variable("x[$i]", VariableOrder, MonomialOrder)
    model.variables[vi] = var
    return vi
end

function _polynomial(variables, func::PolyJuMP.ScalarPolynomialFunction)
    new_variables = VarType[variables[vi] for vi in func.variables]
    return func.polynomial(MP.variables(func.polynomial) => new_variables)
end

function _set(poly::MP.AbstractPolynomialLike, set::MOI.EqualTo)
    return SS.equality(poly, MOI.constant(set))
end

function _set(poly::MP.AbstractPolynomialLike, set::MOI.GreaterThan)
    return SS.PolynomialInequality(poly - MOI.constant(set))
end

function _set(poly::MP.AbstractPolynomialLike, set::MOI.LessThan)
    return SS.PolynomialInequality(MOI.constant(set) - poly)
end

_nineq(::SS.AbstractAlgebraicSet) = 0
_nineq(set) = SS.ninequalities(set)

_num(set, ::Type{<:MOI.EqualTo}) = SS.nequalities(set)
function _num(set, ::Type{<:Union{MOI.LessThan,MOI.GreaterThan}})
    return _nineq(set)
end

function MOI.supports_constraint(
    ::Model{T},
    ::Type{<:PolyJuMP.ScalarPolynomialFunction{T}},
    ::Type{<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T}}},
) where {T}
    return true
end

function MOI.is_valid(
    model::Model,
    ci::MOI.ConstraintIndex{<:PolyJuMP.ScalarPolynomialFunction,S},
) where {S}
    return ci.value in 1:_num(model.set, S)
end

function MOI.add_constraint(
    model::Model{T},
    func::PolyJuMP.ScalarPolynomialFunction{T},
    set::MOI.AbstractScalarSet,
) where {T}
    model.set = model.set ∩ _set(_polynomial(model.variables, func), set)
    i = _num(model.set, typeof(set))
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(i)
end

function MOI.set(
    model::Model,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    model.objective_sense = sense
    return
end

function MOI.supports(
    ::Model{T},
    ::MOI.ObjectiveFunction{<:PolyJuMP.ScalarPolynomialFunction{T}},
) where {T}
    return true
end

function MOI.set(
    model::Model{T},
    ::MOI.ObjectiveFunction,
    func::PolyJuMP.ScalarPolynomialFunction{T},
) where {T}
    model.objective_function = _polynomial(model.variables, func)
    return
end

function _add_to_system(_, lagrangian, ::SS.FullSpace, ::Bool)
    return lagrangian
end

function _add_to_system(
    system,
    lagrangian,
    set::SS.AlgebraicSet,
    maximization::Bool,
)
    n = SS.nequalities(set)
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
    DynamicPolynomials.@polyvar σ[1:PolyJuMP._nineq(set)]
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

function lagrangian_kkt(
    objective_sense::MOI.ObjectiveSense,
    objective_function::MP.AbstractPolynomialLike{T},
    set;
    solver = nothing,
    variables = nothing,
) where {T}
    if isnothing(solver)
        system = SS.AlgebraicSet{T,PolyJuMP.PolyType{T}}()
    else
        I = SS.PolynomialIdeal{T,PolyJuMP.PolyType{T}}()
        system = SS.AlgebraicSet(I, solver)
    end
    if objective_sense == MOI.FEASIBILITY_SENSE
        lagrangian = MA.Zero()
    else
        lagrangian = MA.mutable_copy(objective_function)
    end
    lagrangian = _add_to_system(
        system,
        lagrangian,
        set,
        objective_sense == MOI.MAX_SENSE,
    )
    if !(lagrangian isa MA.Zero)
        ∇x = MP.differentiate(lagrangian, variables)
        for p in ∇x
            SS.add_equality!(system, p)
        end
    end
    return lagrangian, system
end

function lagrangian_kkt(model::Model{T}, solver = nothing) where {T}
    return lagrangian_kkt(
        model.objective_sense,
        model.objective_function,
        model.set;
        solver,
        variables = MP.variables(model),
    )
end
