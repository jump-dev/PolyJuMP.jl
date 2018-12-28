struct ZeroPolynomialSet{BT <: AbstractPolynomialBasis,
                         MT <: AbstractMonomial,
                         MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    basis::Type{BT}
    monomials::MVT
end

struct ZeroPolynomialSetInDomain{BT <: AbstractPolynomialBasis,
                                 DT <: AbstractSemialgebraicSet,
                                 MT <: AbstractMonomial,
                                 MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    basis::Type{BT}
    domain::DT
    monomials::MVT
end
