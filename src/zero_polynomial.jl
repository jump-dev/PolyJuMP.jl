struct ZeroPolynomialSet{BT <: AbstractPolynomialBasis,
                         MT <: AbstractMonomial,
                         MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    basis::Type{BT}
    monomials::MVT
end

struct ZeroPolynomialSetInDomain{DT <: AbstractSemialgebraicSet,
                                 BT <: AbstractPolynomialBasis,
                                 MT <: AbstractMonomial,
                                 MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    domain::DT
    basis::Type{BT}
    monomials::MVT
end
