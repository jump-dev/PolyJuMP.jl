struct ZeroPolynomialSet{DT <: AbstractSemialgebraicSet,
                         BT <: AbstractPolynomialBasis,
                         MT <: AbstractMonomial,
                         MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    domain::DT
    basis::Type{BT}
    monomials::MVT
end
