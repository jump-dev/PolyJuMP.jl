import MathOptInterface as MOI

struct Optimizer{T,O<:MOI.ModelLike} <: MOI.AbstractOptimizer
    model::O
end
function Optimizer(model::MOI.ModelLike)
    return Optimizer{Float64,typeof(model)}(model)
end

function MOI.get(
    model::Optimizer{T},
    attr::MOI.Bridges.ListOfNonstandardBridges,
) where {T}
    list = copy(MOI.get(model.model, attr))
    push!(list, PolyJuMP.Bridges.Constraint.ToPolynomialBridge{T})
    push!(list, PolyJuMP.Bridges.Objective.ToPolynomialBridge{T})
    return list
end

MOI.is_empty(model::Optimizer) = MOI.is_empty(model.model)
MOI.empty!(model::Optimizer) = MOI.empty!(model.model)

MOI.add_variable(model::Optimizer) = MOI.add_variable(model.model)

function MOI.supports_add_constrained_variable(
    model::Optimizer,
    ::Type{S},
) where {S<:MOI.AbstractScalarSet}
    return MOI.supports_add_constrained_variable(model.model, S)
end

function MOI.supports_add_constrained_variables(
    model::Optimizer,
    ::Type{MOI.Reals},
)
    return MOI.supports_add_constrained_variables(model.model, MOI.Reals)
end

function MOI.supports_add_constrained_variables(
    model::Optimizer,
    ::Type{S},
) where {S<:MOI.AbstractVectorSet}
    return MOI.supports_add_constrained_variables(model.model, S)
end

function MOI.supports(model::Optimizer, attr::MOI.AbstractModelAttribute)
    return MOI.supports(model.model, attr)
end

function MOI.set(model::Optimizer, attr::MOI.AbstractModelAttribute, value)
    MOI.set(model.model, attr, value)
end

function MOI.get(model::Optimizer, attr::MOI.AbstractModelAttribute)
    return MOI.get(model.model, attr)
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
    return MOI.add_constraint(model.model, func, set)
end

function MOI.supports_incremental_interface(model::Optimizer)
    return MOI.supports_incremental_interface(model.model)
end
function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

function MOI.get(model::Optimizer, attr::MOI.SolverName)
    name = MOI.get(model.model, attr)
    return "PolyJuMP.QCQP with $name"
end
