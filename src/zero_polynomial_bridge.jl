struct ZeroPolynomialBridge{T, F <: MOI.AbstractVectorFunction,
                            MT <: AbstractMonomial,
                            MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    zero_constraint::MOI.ConstraintIndex{F, MOI.Zeros}
    monomials::MVT
end

function ZeroPolynomialBridge{T, F}(model::MOI.ModelLike,
                                    f::MOI.AbstractVectorFunction,
                                    s::ZeroPolynomialSet{<:MonomialBasis}) where {T, F}
    @assert MOI.output_dimension(f) == s.monomials
    zero_constraint = MOI.add_constraint(model, f,
                                         MOI.Zeros(length(s.monomials)))
    return ZeroPolynomialBridge{T, F}(zero_constraint, s.monomials)
end

function MOI.supports_constraint(::Type{ZeroPolynomialBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:ZeroPolynomialSet}) where T
    return true
end
function added_constraint_types(::Type{<:ZeroPolynomialBridge{T, F}}) where {T, F}
    return [(F, MOI.Zeros)]
end
function concrete_bridge_type(::Type{<:ZeroPolynomialBridge{T}},
                              F::Type{<:MOI.AbstractVectorFunction},
                              ::Type{ZeroPolynomialSet{<:MonomialBasis, MT, MVT}}) where {T, MT, MVT}
    return ZeroPolynomialBridge{T, F, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(::ZeroPolynomialBridge{T, F},
                 ::MOI.NumberOfConstraints{F, MOI.Zeros}) where {T, F}
    return 1
end
function MOI.get(b::ZeroPolynomialBridge{T, F},
                 ::MOI.ListOfConstraintIndices{F, MOI.Zeros}) where {T, F}
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, c::ZeroPolynomialBridge)
    MOI.delete(model, c.zero_constraint)
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 attr::Union{MOI.ConstraintPrimal, MOI.ConstraintDual},
                 c::ZeroPolynomialBridge)
    return MOI.get(model, attr, c.zero_constraint)
end
