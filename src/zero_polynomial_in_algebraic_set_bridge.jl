struct ZeroPolynomialInAlgebraicSetBridge{T, F <: MOI.AbstractVectorFunction,
                                          BT <: AbstractPolynomialBasis,
                                          MT <: AbstractMonomial,
                                          MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    zero_constraint::MOI.ConstraintIndex{F, ZeroPolynomialSet{BT, MT, MVT}}
end

function ZeroPolynomialInAlgebraicSetBridge{T, F}(model::MOI.ModelLike,
                                                  f::MOI.AbstractVectorFunction,
                                                  s::ZeroPolynomialSetInDomain{<:AbstractAlgebraicSet}) where {T, F}
    p = polynomial(MOI.Utilities.eachscalar(f), s.monomials)
    r = rem(p, ideal(s.domain))
    zero_constraint = MOI.add_constraint(model, coefficients(r),
                                         ZeroPolynomialSet(s.basis,
                                                           monomials(r)))
    return ZeroPolynomialInAlgebraicSetBridge{T, F}(zero_constraint)
end


function MOI.supports_constraint(::Type{ZeroPolynomialInAlgebraicSetBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:ZeroPolynomialSetInDomain{<:AbstractAlgebraicSet}}) where T
    return true
end
function added_constraint_types(::Type{<:ZeroPolynomialInAlgebraicSetBridge{T, F, BT, MT, MVT}}) where {T, F, BT, MT, MVT}
    return [(F, ZeroPolynomialSet{BT, MT, MVT})]
end
function concrete_bridge_type(::Type{<:ZeroPolynomialInAlgebraicSetBridge{T}},
                              F::Type{<:MOI.AbstractVectorFunction},
                              ::Type{ZeroPolynomialSet{BT, MT, MVT}}) where {T, BT, MT, MVT}
    G = MOI.Utilities.promote_operation(-, T, F, F)
    return ZeroPolynomialInAlgebraicSetBridge{T, G, BT, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(::ZeroPolynomialInAlgebraicSetBridge{T, F, BT, MT, MVT},
                 ::MOI.NumberOfConstraints{F, ZeroPolynomialSet{BT, MT, MVT}}) where {T, F, BT, MT, MVT}
    return 1
end
function MOI.get(::ZeroPolynomialInAlgebraicSetBridge{T, F, BT, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, ZeroPolynomialSet{BT, MT, MVT}}) where {T, F, BT, MT, MVT}
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, c::ZeroPolynomialInAlgebraicSetBridge)
    MOI.delete(model, c.zero_constraint)
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 attr::Union{MOI.ConstraintPrimal, MOI.ConstraintDual},
                 c::ZeroPolynomialInAlgebraicSetBridge)
    return MOI.get(model, attr, c.zero_constraint)
end
