struct ZeroPolynomialInAlgebraicSetBridge{T, F <: MOI.AbstractVectorFunction,
                                          BT <: AbstractPolynomialBasis,
                                          MT <: AbstractMonomial,
                                          MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    zero_constraint::MOI.ConstraintIndex{F, ZeroPolynomialSet{FullSpace, BT, MT, MVT}}
end

function ZeroPolynomialInAlgebraicSetBridge{T, F, BT, MT, MVT}(model::MOI.ModelLike,
                                                               f::MOI.AbstractVectorFunction,
                                                               s::ZeroPolynomialSet{<:AbstractAlgebraicSet}) where {T, F, BT, MT, MVT}
    p = polynomial(collect(MOI.Utilities.eachscalar(f)), s.monomials)
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    r = convert(typeof(p), rem(p, ideal(s.domain)))
    zero_constraint = MOI.add_constraint(model, MOIU.vectorize(coefficients(r)),
                                         ZeroPolynomialSet(FullSpace(), s.basis,
                                                           monomials(r)))
    return ZeroPolynomialInAlgebraicSetBridge{T, F, BT, MT, MVT}(zero_constraint)
end


function MOI.supports_constraint(::Type{ZeroPolynomialInAlgebraicSetBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:ZeroPolynomialSet{<:AbstractAlgebraicSet}}) where T
    return true
end
function MOIB.added_constraint_types(::Type{<:ZeroPolynomialInAlgebraicSetBridge{T, F, BT, MT, MVT}}) where {T, F, BT, MT, MVT}
    return [(F, ZeroPolynomialSet{FullSpace, BT, MT, MVT})]
end
function MOIB.concrete_bridge_type(::Type{<:ZeroPolynomialInAlgebraicSetBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{<:ZeroPolynomialSet{<:AbstractAlgebraicSet,
                                                              BT, MT, MVT}}) where {T, BT, MT, MVT}
    G = MOI.Utilities.promote_operation(-, T, F, F)
    return ZeroPolynomialInAlgebraicSetBridge{T, G, BT, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(::ZeroPolynomialInAlgebraicSetBridge{T, F, BT, MT, MVT},
                 ::MOI.NumberOfConstraints{F, ZeroPolynomialSet{FullSpace, BT, MT, MVT}}) where {T, F, BT, MT, MVT}
    return 1
end
function MOI.get(::ZeroPolynomialInAlgebraicSetBridge{T, F, BT, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, ZeroPolynomialSet{FullSpace, BT, MT, MVT}}) where {T, F, BT, MT, MVT}
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
