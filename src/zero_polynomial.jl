struct ZeroPolynomialSet{DT <: AbstractSemialgebraicSet,
                         BT <: MB.AbstractPolynomialBasis,
                         MT <: AbstractMonomial,
                         MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    domain::DT
    basis::Type{BT}
    monomials::MVT
end

# `x`-in-`PlusMinusSet(set)` iff `x`-in-`set` and `-x`-in-`set`.
struct PlusMinusSet{ST <: MOI.AbstractVectorSet} <: MOI.AbstractVectorSet
    set::ST
end
