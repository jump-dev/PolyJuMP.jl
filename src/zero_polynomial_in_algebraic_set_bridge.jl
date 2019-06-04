using LinearAlgebra

struct ZeroPolynomialInAlgebraicSetBridge{T, F <: MOI.AbstractVectorFunction,
                                          BT <: AbstractPolynomialBasis,
                                          DT <: AbstractSemialgebraicSet,
                                          MT <: AbstractMonomial,
                                          MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    zero_constraint::MOI.ConstraintIndex{F, ZeroPolynomialSet{FullSpace, BT, MT, MVT}}
    domain::DT
    monomials::MVT
end

function ZeroPolynomialInAlgebraicSetBridge{T, F, BT, DT, MT, MVT}(model::MOI.ModelLike,
                                                                   f::MOI.AbstractVectorFunction,
                                                                   s::ZeroPolynomialSet{<:AbstractAlgebraicSet}) where {T, F <: MOI.AbstractVectorFunction,
                                                                                                                        BT <: AbstractPolynomialBasis,
                                                                                                                        DT <: AbstractSemialgebraicSet,
                                                                                                                        MT <: AbstractMonomial,
                                                                                                                        MVT <: AbstractVector{MT}}
    p = polynomial(collect(MOI.Utilities.eachscalar(f)), s.monomials)
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    r = convert(typeof(p), rem(p, ideal(s.domain)))
    zero_constraint = MOI.add_constraint(model, MOIU.vectorize(coefficients(r)),
                                         ZeroPolynomialSet(FullSpace(), s.basis,
                                                           monomials(r)))
    return ZeroPolynomialInAlgebraicSetBridge{T, F, BT, DT, MT, MVT}(zero_constraint, s.domain, s.monomials)
end


function MOI.supports_constraint(::Type{ZeroPolynomialInAlgebraicSetBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:ZeroPolynomialSet{<:AbstractAlgebraicSet}}) where T
    return true
end
function MOIB.added_constraint_types(::Type{<:ZeroPolynomialInAlgebraicSetBridge{T, F, BT, DT, MT, MVT}}) where {T, F, BT, DT, MT, MVT}
    return [(F, ZeroPolynomialSet{FullSpace, BT, MT, MVT})]
end
function MOIB.concrete_bridge_type(::Type{<:ZeroPolynomialInAlgebraicSetBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{<:ZeroPolynomialSet{DT, BT, MT, MVT}}) where {T, BT, DT<:AbstractAlgebraicSet, MT, MVT}
    G = MOI.Utilities.promote_operation(-, T, F, F)
    return ZeroPolynomialInAlgebraicSetBridge{T, G, BT, DT, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(::ZeroPolynomialInAlgebraicSetBridge{T, F, BT, DT, MT, MVT},
                 ::MOI.NumberOfConstraints{F, ZeroPolynomialSet{FullSpace, BT, MT, MVT}}) where {T, F, BT, DT, MT, MVT}
    return 1
end
function MOI.get(b::ZeroPolynomialInAlgebraicSetBridge{T, F, BT, DT, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, ZeroPolynomialSet{FullSpace, BT, MT, MVT}}) where {T, F, BT, DT, MT, MVT}
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, c::ZeroPolynomialInAlgebraicSetBridge)
    MOI.delete(model, c.zero_constraint)
end

# Attributes, Bridge acting as a constraint

# TODO ConstraintPrimal

# Let A be the linear map corresponding to A(p) = rem(p, ideal(set.domain))
# We want to compute A*(μ) (where A* is the conjugate of A); see
# slide 9 of https://www.youtube.com/watch?v=C8dHxJCUHYw.
# That is, for every monomial of the basis `mono`, we want to compute
# `⟨mono, A*(μ)⟩`. By definition of the conjugacy, we have `⟨mono, A*(μ)⟩ = ⟨A(mono), μ⟩`
# so we can compute compute it by `⟨rem(mono, ideal(set.domain)), μ⟩`
function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::ZeroPolynomialInAlgebraicSetBridge)
    dual = MOI.get(model, attr, bridge.zero_constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.zero_constraint)
    μ = measure(dual, set.monomials)
    I = ideal(bridge.domain)
    return [dot(rem(mono, I), μ) for mono in bridge.monomials]
end
function MOI.get(model::MOI.ModelLike, attr::MomentsAttribute,
                 bridge::ZeroPolynomialInAlgebraicSetBridge)
    return MOI.get(model, attr, bridge.zero_constraint)
end
