struct ZeroPolynomialInSemialgebraicSetBridge{T, F <: MOI.AbstractVectorFunction,
                                          BT <: AbstractPolynomialBasis,
                                          MT <: AbstractMonomial,
                                          MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    zero_constraint::MOI.ConstraintIndex{F, ZeroPolynomialSet{FullSpace, BT, MT, MVT}}
end

function ZeroPolynomialInSemialgebraicSetBridge{T, F}(model::MOI.ModelLike,
                                                  f::MOI.AbstractVectorFunction,
                                                  s::ZeroPolynomialSet{<:BasicSemialgebraicSet}) where {T, F}
    p = polynomial(MOI.Utilities.eachscalar(f), s.monomials)
    r = rem(p, ideal(s.domain))
    zero_constraint = MOI.add_constraint(model, coefficients(r),
                                         ZeroPolynomialSet(FullSpace(), s.basis,
                                                           monomials(r)))
    return ZeroPolynomialInSemialgebraicSetBridge{T, F}(zero_constraint)
end


function MOI.supports_constraint(::Type{ZeroPolynomialInSemialgebraicSetBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:ZeroPolynomialSet{<:BasicSemialgebraicSet}}) where T
    return true
end
function added_constraint_types(::Type{<:ZeroPolynomialInSemialgebraicSetBridge{T, F, BT, MT, MVT}}) where {T, F, BT, MT, MVT}
    return [(F, ZeroPolynomialSet{FullSpace, BT, MT, MVT})]
end
function concrete_bridge_type(::Type{<:ZeroPolynomialInSemialgebraicSetBridge{T}},
                              F::Type{<:MOI.AbstractVectorFunction},
                              ::Type{ZeroPolynomialSet{<:BasicSemialgebraicSet,
                                                       BT, MT, MVT}}) where {T, BT, MT, MVT}
    G = MOI.Utilities.promote_operation(-, T, F, F)
    return ZeroPolynomialInSemialgebraicSetBridge{T, G, BT, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(::ZeroPolynomialInSemialgebraicSetBridge{T, F, BT, MT, MVT},
                 ::MOI.NumberOfConstraints{F, ZeroPolynomialSet{FullSpace, BT, MT, MVT}}) where {T, F, BT, MT, MVT}
    return 1
end
function MOI.get(::ZeroPolynomialInSemialgebraicSetBridge{T, F, BT, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, ZeroPolynomialSet{FullSpace, BT, MT, MVT}}) where {T, F, BT, MT, MVT}
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, c::ZeroPolynomialInSemialgebraicSetBridge)
    MOI.delete(model, c.zero_constraint)
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 attr::Union{MOI.ConstraintPrimal, MOI.ConstraintDual},
                 c::ZeroPolynomialInSemialgebraicSetBridge)
    return MOI.get(model, attr, c.zero_constraint)
end
