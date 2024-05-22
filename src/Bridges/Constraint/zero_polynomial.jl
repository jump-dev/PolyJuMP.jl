struct ZeroPolynomialBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    B<:SA.ExplicitBasis,
} <: MOI.Bridges.Constraint.AbstractBridge
    zero_constraint::MOI.ConstraintIndex{F,MOI.Zeros}
    basis::B
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ZeroPolynomialBridge{T,F,B}},
    model::MOI.ModelLike,
    f::MOI.AbstractVectorFunction,
    s::PolyJuMP.ZeroPolynomialSet{SS.FullSpace,<:MB.FullBasis{MB.Monomial},B},
) where {T,F,B}
    @assert MOI.output_dimension(f) == length(s.basis)
    zero_constraint =
        MOI.add_constraint(model, f, MOI.Zeros(length(s.basis)))
    return ZeroPolynomialBridge{T,F,typeof(s.basis)}(zero_constraint, s.basis)
end

function MOI.supports_constraint(
    ::Type{<:ZeroPolynomialBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:PolyJuMP.ZeroPolynomialSet{SS.FullSpace}},
) where {T}
    return true
end
function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:ZeroPolynomialBridge},
)
    return Tuple{Type}[]
end
function MOI.Bridges.added_constraint_types(
    ::Type{<:ZeroPolynomialBridge{T,F}},
) where {T,F}
    return [(F, MOI.Zeros)]
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ZeroPolynomialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{
        <:PolyJuMP.ZeroPolynomialSet{SS.FullSpace,<:MB.FullBasis{MB.Monomial},B},
    },
) where {T,B<:MB.SubBasis{MB.Monomial}}
    return ZeroPolynomialBridge{T,F,B}
end

# Attributes, Bridge acting as an model
function MOI.get(
    ::ZeroPolynomialBridge{T,F},
    ::MOI.NumberOfConstraints{F,MOI.Zeros},
) where {T,F}
    return 1
end
function MOI.get(
    b::ZeroPolynomialBridge{T,F},
    ::MOI.ListOfConstraintIndices{F,MOI.Zeros},
) where {T,F}
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::ZeroPolynomialBridge)
    return MOI.delete(model, bridge.zero_constraint)
end

# Attributes, Bridge acting as a constraint
function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::ZeroPolynomialBridge{T,F,MB.SubBasis{MB.Monomial,M,V}},
) where {T,F,M,V}
    return PolyJuMP.ZeroPolynomialSet(
        SS.FullSpace(),
        MB.FullBasis{MB.Monomial,M}(),
        bridge.basis,
    )
end
function MOI.get(
    model::MOI.ModelLike,
    attr::Union{MOI.ConstraintPrimal,MOI.ConstraintDual},
    bridge::ZeroPolynomialBridge,
)
    return MOI.get(model, attr, bridge.zero_constraint)
end
function MOI.get(
    model::MOI.ModelLike,
    attr::PolyJuMP.MomentsAttribute,
    bridge::ZeroPolynomialBridge,
)
    values = MOI.get(model, MOI.ConstraintDual(attr.result_index), bridge)
    return MM.measure(values, bridge.basis)
end
