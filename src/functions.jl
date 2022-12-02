"""
    struct ScalarPolynomialFunction{P<:MP.AbstractPolynomialLike} <: MOI.AbstractScalarFunction
        p::P
        variables::Vector{MOI.VariableIndex}
    end

Defines the polynomial function of the variables `variables` where the variable
`variables(p)[i]` corresponds to `variables[i]`.
"""
struct ScalarPolynomialFunction{T,P<:AbstractPolynomial{T}} <: MOI.AbstractScalarFunction
    polynomial::P
    variables::Vector{MOI.VariableIndex}
end
