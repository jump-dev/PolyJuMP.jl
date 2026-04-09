struct ZeroPolynomialSet{
    DT<:SS.AbstractSemialgebraicSet,
    Z<:SA.AbstractBasis,
    B<:SA.ExplicitBasis,
} <: MOI.AbstractVectorSet
    domain::DT
    zero_basis::Z
    basis::B
end

MOI.dimension(set::ZeroPolynomialSet) = length(set.basis)
Base.copy(set::ZeroPolynomialSet) = set

# `x`-in-`PlusMinusSet(set)` iff `x`-in-`set` and `-x`-in-`set`.
struct PlusMinusSet{ST<:MOI.AbstractVectorSet} <: MOI.AbstractVectorSet
    set::ST
end
MOI.dimension(set::PlusMinusSet) = MOI.dimension(set.set)
Base.copy(set::PlusMinusSet) = set
