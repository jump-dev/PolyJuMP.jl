struct PlusMinusBridge{T, F <: MOI.AbstractVectorFunction,
                       G <: MOI.AbstractVectorFunction,
                       ST <: MOI.AbstractVectorSet} <: MOIB.AbstractBridge
    plus::MOI.ConstraintIndex{F, ST}
    minus::MOI.ConstraintIndex{G, ST}
end

function PlusMinusBridge{T, F, G, ST}(model::MOI.ModelLike,
                                      f::MOI.AbstractVectorFunction,
                                      set::PlusMinusSet{ST}) where {T, F, G, ST}
    plus = MOI.add_constraint(model, f, set)
    minus = MOI.add_constraint(model, operate(-, T, f), set)
    return PlusMinusBridge{T, F, G, ST}(plus, minus)
end


function MOI.supports_constraint(::Type{PlusMinusBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:PlusMinusSet}) where T
    return true
end
function MOIB.added_constraint_types(::Type{<:PlusMinusBridge{T, F, F, ST}}) where {T, F, ST}
    return [(F, ST)]
end
function MOIB.added_constraint_types(::Type{<:PlusMinusBridge{T, F, G, ST}}) where {T, F, G, ST}
    return [(F, ST), (G, ST)]
end
function MOIB.concrete_bridge_type(::Type{<:PlusMinusBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{<:PlusMinusSet{ST}}) where {T, ST}
    G = MOI.Utilities.promote_operation(-, T, F)
    return PlusMinusBridge{T, F, G, ST}
end

# Attributes, Bridge acting as an model
function MOI.get(::PlusMinusBridge{T, F, F, ST},
                 ::MOI.NumberOfConstraints{F, ST}) where {T, F, ST}
    return 1
end
function MOI.get(::PlusMinusBridge{T, F, G, ST},
                 ::MOI.NumberOfConstraints{F, ST}) where {T, F, G, ST}
    return 1
end
function MOI.get(::PlusMinusBridge{T, F, G, ST},
                 ::MOI.NumberOfConstraints{G, ST}) where {T, F, G, ST}
    return 1
end
function MOI.get(bridge::PlusMinusBridge{T, F, F, ST},
                 ::MOI.ListOfConstraintIndices{F, ST}) where {T, F, ST}
    return [bridge.plus, bridge.minus]
end
function MOI.get(bridge::PlusMinusBridge{T, F, G, ST},
                 ::MOI.ListOfConstraintIndices{F, ST}) where {T, F, G, ST}
    return [bridge.plus]
end
function MOI.get(bridge::PlusMinusBridge{T, F, G, ST},
                 ::MOI.ListOfConstraintIndices{G, ST}) where {T, F, G, ST}
    return [bridge.minus]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::PlusMinusBridge)
    MOI.delete(model, bridge.plus)
    MOI.delete(model, bridge.minus)
end
