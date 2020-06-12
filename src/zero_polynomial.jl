struct ZeroPolynomialSet{DT <: AbstractSemialgebraicSet,
                         BT <: MB.AbstractPolynomialBasis
                        } <: MOI.AbstractVectorSet
    domain::DT
    basis::BT
end

# `x`-in-`PlusMinusSet(set)` iff `x`-in-`set` and `-x`-in-`set`.
struct PlusMinusSet{ST <: MOI.AbstractVectorSet} <: MOI.AbstractVectorSet
    set::ST
end
