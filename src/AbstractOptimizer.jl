mutable struct AbstractOptimizer{T} <: MOI.AbstractOptimizer end

function invalidate_solution! end

function MOI.get(
    ::AbstractOptimizer{T},
    ::MOI.Bridges.ListOfNonstandardBridges{T},
) where {T}
    return [
        PolyJuMP.Bridges.Constraint.ToPolynomialBridge{T},
        PolyJuMP.Bridges.Objective.ToPolynomialBridge{T},
    ]
end

function MOI.is_empty(model::AbstractOptimizer)
    return isempty(model.variables) &&
           model.objective_sense == MOI.FEASIBILITY_SENSE &&
           isnothing(model.objective_function) &&
           model.set isa SemialgebraicSets.FullSpace
end

function MOI.empty!(model::AbstractOptimizer)
    empty!(model.variables)
    model.objective_sense = MOI.FEASIBILITY_SENSE
    model.objective_function = nothing
    model.set = SemialgebraicSets.FullSpace()
    invalidate_solution!(model)
    return
end

function MOI.is_valid(model::AbstractOptimizer, vi::MOI.VariableIndex)
    return in(vi.value, 1:length(model.variables))
end

function MOI.add_variable(model::AbstractOptimizer)
    i = length(model.variables) + 1
    vi = MOI.VariableIndex(i)
    var = DynamicPolynomials.Variable("x[$i]", VariableOrder, MonomialOrder)
    model.variables[vi] = var
    invalidate_solution!(model)
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
    ::AbstractOptimizer{T},
    ::Type{<:PolyJuMP.ScalarPolynomialFunction{T}},
    ::Type{<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T}}},
) where {T}
    return true
end

function MOI.is_valid(
    model::AbstractOptimizer,
    ci::MOI.ConstraintIndex{<:PolyJuMP.ScalarPolynomialFunction,S},
) where {S}
    return ci.value in 1:_num(model.set, S)
end

function MOI.add_constraint(
    model::AbstractOptimizer{T},
    func::PolyJuMP.ScalarPolynomialFunction{T},
    set::MOI.AbstractScalarSet,
) where {T}
    model.set = model.set ∩ _set(_polynomial(model.variables, func), set)
    i = _num(model.set, typeof(set))
    invalidate_solution!(model)
    return MOI.ConstraintIndex{typeof(func),typeof(set)}(i)
end

function MOI.set(
    model::AbstractOptimizer,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    model.objective_sense = sense
    invalidate_solution!(model)
    return
end

function MOI.supports(
    ::AbstractOptimizer{T},
    ::MOI.ObjectiveFunction{<:PolyJuMP.ScalarPolynomialFunction{T}},
) where {T}
    return true
end
function MOI.set(
    model::AbstractOptimizer{T},
    ::MOI.ObjectiveFunction,
    func::PolyJuMP.ScalarPolynomialFunction{T},
) where {T}
    model.objective_function = _polynomial(model.variables, func)
    invalidate_solution!(model)
    return
end

MOI.supports_incremental_interface(::AbstractOptimizer) = true
function MOI.copy_to(dest::AbstractOptimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end
