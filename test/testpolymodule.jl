module TestPolyModule

using LinearAlgebra
using MathOptInterface
const MOI = MathOptInterface
using JuMP
using PolyJuMP
using MultivariatePolynomials
using SemialgebraicSets

struct NonNeg{BT <: PolyJuMP.AbstractPolynomialBasis,
              DT <: SemialgebraicSets.AbstractSemialgebraicSet,
              MT <: MultivariatePolynomials.AbstractMonomial,
              MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    basis::Type{BT}
    domain::DT
    monomials::MVT
    kwargs
end
struct TestNonNeg <: PolyJuMP.PolynomialSet end

function JuMP.moi_set(cone::TestNonNeg,
                      monos::AbstractVector{<:AbstractMonomial};
                      domain::AbstractSemialgebraicSet=FullSpace(),
                      basis=MonomialBasis, kwargs...)
    return NonNeg(basis, domain, monos, kwargs)
end


function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::TestNonNeg; kwargs...)
    monos = monomials(p)
    set = JuMP.moi_set(s, monos; kwargs...)
    return JuMP.build_constraint(_error, coefficients(p), set)
end

struct PosDefMatrix{BT <: PolyJuMP.AbstractPolynomialBasis,
                    DT <: SemialgebraicSets.AbstractSemialgebraicSet,
                    VT <: MultivariatePolynomials.AbstractVariable,
                    MT <: MultivariatePolynomials.AbstractMonomial,
                    MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    basis::Type{BT}
    domain::DT
    y::Vector{VT}
    monomials::MVT
    kwargs
end
struct TestPosDefMatrix <: PolyJuMP.PolynomialSet end
function JuMP.build_constraint(_error::Function,
                               p::Matrix{<:AbstractPolynomialLike},
                               s::TestPosDefMatrix;
                               basis=PolyJuMP.MonomialBasis,
                               domain=FullSpace(), kwargs...)
    n = LinearAlgebra.checksquare(p)
    y = [similarvariable(eltype(p), gensym()) for i in 1:n]
    q = dot(y, p * y)
    set = PosDefMatrix(basis, domain, y, monomials(q), kwargs)
    return JuMP.build_constraint(_error, coefficients(q), set)
end


function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, TestNonNeg)
    PolyJuMP.setdefault!(data, PolyJuMP.PosDefPolyMatrix, TestPosDefMatrix)
end

end
