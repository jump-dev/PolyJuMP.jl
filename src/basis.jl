export FixedPolynomialBasis, ScaledMonomialBasis, MonomialBasis

"""
    abstract type AbstractPolynomialBasis end

Polynomial basis of a subspace of the polynomials [Section~3.1.5, BPT12].

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.
"""
abstract type AbstractPolynomialBasis end

"""
    struct FixedPolynomialBasis{PT<:MultivariatePolynomials.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
        polynomials::PV
    end

Polynomial basis with the polynomials of the vector `polynomials`.
For instance, `FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])` is the Chebyshev
polynomial basis for cubic polynomials in the variable `x`.
"""
struct FixedPolynomialBasis{PT<:MultivariatePolynomials.AbstractPolynomialLike, PV<:AbstractVector{PT}} <: AbstractPolynomialBasis
    polynomials::PV
end

function MultivariatePolynomials.polynomialtype(mb::FixedPolynomialBasis{PT}, T::Type) where PT
    C = MultivariatePolynomials.coefficienttype(PT)
    U = typeof(zero(C) * zero(T) + zero(C) * zero(T))
    MultivariatePolynomials.polynomialtype(PT, U)
end
function MultivariatePolynomials.polynomial(f::Function, fpb::FixedPolynomialBasis)
    sum(ip -> f(ip[1]) * ip[2], enumerate(fpb.polynomials))
end

"""
    struct MonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
        monomials::MV
    end

Monomial basis with the monomials of the vector `monomials`.
For instance, `MonomialBasis([1, x, y, x^2, x*y, y^2])` is the monomial basis
for the subspace of quadratic polynomials in the variables `x`, `y`.
"""
struct MonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
    monomials::MV
end
MonomialBasis(monomials) = MonomialBasis(monovec(monomials))

MultivariatePolynomials.polynomialtype(mb::MonomialBasis{MT}, T::Type) where MT = MultivariatePolynomials.polynomialtype(MT, T)
MultivariatePolynomials.polynomial(f::Function, mb::MonomialBasis) = polynomial(f, mb.monomials)
function MultivariatePolynomials.coefficients(p, ::Type{<:MonomialBasis})
    return coefficients(p)
end

"""
    struct ScaledMonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
        monomials::MV
    end

*Scaled monomial basis* (see [Section 3.1.5, BPT12]) with the monomials of the vector `monomials`.
Given a monomial ``x^\\alpha = x_1^{\\alpha_1} \\cdots x_n^{\\alpha_n}`` of degree ``d = \\sum_{i=1}^n \\alpha_i``,
the corresponding polynomial of the basis is
```math
{d \\choose \\alpha}^{\\frac{1}{2}} x^{\\alpha} \\quad \\text{ where } \\quad
{d \\choose \\alpha} = \\frac{d!}{\\alpha_1! \\alpha_2! \\cdots \\alpha_n!}.
```

For instance, create a polynomial with the basis ``[xy^2, xy]`` creates the polynomial
``\\sqrt{3} a xy^2 + \\sqrt{2} b xy`` where `a` and `b` are new JuMP decision variables.
Constraining the polynomial ``axy^2 + bxy`` to be zero with the scaled monomial basis constrains
`a/√3` and `b/√2` to be zero.

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.
"""
struct ScaledMonomialBasis{MT<:MultivariatePolynomials.AbstractMonomial, MV<:AbstractVector{MT}} <: AbstractPolynomialBasis
    monomials::MV
end
ScaledMonomialBasis(monomials) = ScaledMonomialBasis(monovec(monomials))

MultivariatePolynomials.polynomialtype(mb::ScaledMonomialBasis{MT}, T::Type) where MT = MultivariatePolynomials.polynomialtype(MT, promote_type(T, Float64))
scaling(m::MultivariatePolynomials.AbstractMonomial) = √(factorial(degree(m)) / prod(factorial, exponents(m)))
MultivariatePolynomials.polynomial(f::Function, mb::ScaledMonomialBasis) = polynomial(i -> scaling(mb.monomials[i]) * f(i), mb.monomials)
unscale_coef(t::MultivariatePolynomials.AbstractTerm) = coefficient(t) / scaling(monomial(t))
function MultivariatePolynomials.coefficients(p, ::Type{<:ScaledMonomialBasis})
    return unscale_coef.(terms(p))
end
