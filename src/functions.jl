function Base.:*(f::MOI.Utilities.ScalarLike, p::AbstractPolynomialLike)
    return MultivariatePolynomials.multconstant(f, p)
end
function Base.:*(p::AbstractPolynomialLike, f::MOI.Utilities.ScalarLike)
    return MultivariatePolynomials.multconstant(p, f)
end
