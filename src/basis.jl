"""
    abstract type AbstractPolynomialBasis end

Polynomial basis of a subspace of the polynomials [Section~3.1.5, BPT12].

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, 2012.
"""
abstract type AbstractPolynomialBasis end

"""
    struct FixedPolynomialBasis{PT<:MultivariatePolynomials.AbstractPolynomial, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
        p::PV
    end

Polynomial basis with the polynomials of the vector `p`.
For instance, `FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])` is the Chebyshev
polynomial basis for cubic polynomials in the variable `x`.
"""
struct FixedPolynomialBasis{PT<:MultivariatePolynomials.AbstractPolynomial, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
    p::PV
end

function MultivariatePolynomials.polynomialtype(mb::FixedPolynomialBasis{PT}, T::Type) where PT
    C = MultivariatePolynomials.coefficienttype(PT)
    U = typeof(zero(C) * zero(T) + zero(C) * zero(T))
    MultivariatePolynomials.polynomialtype(PT, U)
end

"""
    struct MonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
        x::MV
    end

Monomial basis with the monomials of the vector `x`.
For instance, `MonomialBasis([1, x, y, x^2, x*y, y^2])` is the monomial basis
for the subspace of quadratic polynomials in the variables `x`, `y`.
"""
struct MonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
    x::MV
end
MonomialBasis(x::AbstractVector{MT}) where {MT<:MultivariatePolynomials.AbstractMonomial} = MonomialBasis{MT, typeof(x)}(x)
MonomialBasis(x) = MonomialBasis(monovec(x))

MultivariatePolynomials.polynomialtype(mb::MonomialBasis{MT}, T::Type) where MT = MultivariatePolynomials.polynomialtype(MT, T)
MultivariatePolynomials.polynomial(f::Function, mb::MonomialBasis) = polynomial(f, mb.x)

abstract type AbstractMeasureBasis end
