# TODO remove this type piracy when we require MOI v0.9 as we need this fix: https://github.com/jump-dev/MathOptInterface.jl/pull/800
function Base.:*(f::MOI.Utilities.ScalarLike, p::AbstractPolynomialLike)
    return MultivariatePolynomials.multconstant(f, p)
end
function Base.:*(p::AbstractPolynomialLike, f::MOI.Utilities.ScalarLike)
    return MultivariatePolynomials.multconstant(p, f)
end
