struct ZeroPolynomialSet{
    DT<:SS.AbstractSemialgebraicSet,
    BT<:MB.AbstractPolynomialBasis,
    MT<:MP.AbstractMonomial,
    MVT<:AbstractVector{MT},
} <: MOI.AbstractVectorSet
    domain::DT
    basis::Type{BT}
    monomials::MVT
    function ZeroPolynomialSet(
        domain::SS.AbstractSemialgebraicSet,
        basis::Type{BT},
        monomials::AbstractVector{MT},
    ) where {BT<:MB.AbstractPolynomialBasis,MT<:MP.AbstractMonomial}
        # For terms, `monomials` is `OneOrZeroElementVector`
        # so we convert it with `monomial_vector`
        # Later, we'll use `MP.MonomialBasis` which is going to do that anyway
        vec = _lazy_monomial_vector(monomials)
        return new{typeof(domain),BT,MT,typeof(vec)}(domain, basis, vec)
    end
end

function _lazy_monomial_vector(monomials)
    MVT = MP.monomial_vector_type(typeof(monomials))
    if MVT != typeof(monomials)
        vec = MP.monomial_vector(monomials)
        @assert typeof(vec) == MVT
        return vec
    end
    return monomials
end

MOI.dimension(set::ZeroPolynomialSet) = length(set.basis)
Base.copy(set::ZeroPolynomialSet) = set

# `x`-in-`PlusMinusSet(set)` iff `x`-in-`set` and `-x`-in-`set`.
struct PlusMinusSet{ST<:MOI.AbstractVectorSet} <: MOI.AbstractVectorSet
    set::ST
end
MOI.dimension(set::PlusMinusSet) = MOI.dimension(set.set)
Base.copy(set::PlusMinusSet) = set
