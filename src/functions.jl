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

function Base.copy(func::ScalarPolynomialFunction)
    return ScalarPolynomialFunction(MA.copy_if_mutable(func.polynomial), copy(func.variables))
end

function MOI.Utilities.canonicalize!(::ScalarPolynomialFunction) end

function MOI.Utilities.is_coefficient_type(
    ::Type{<:ScalarPolynomialFunction{T}},
    ::Type{S},
) where {S,T}
    return S === T
end

function MOI.Utilities.promote_operation(
    ::typeof(-),
    ::Type{T},
    F::Type{ScalarPolynomialFunction{T,P}},
) where {T,P}
    return F
end

function MOI.Utilities.promote_operation(
    ::typeof(-),
    ::Type{T},
    F::Type{ScalarPolynomialFunction{T,P}},
    ::Type{<:Union{T,MOI.Utilities.ScalarLike{T}}},
) where {T,P}
    return F
end

# FIXME
function MOI.Utilities.promote_operation(
    ::typeof(vcat),
    ::Type{T},
    ::Type{ScalarPolynomialFunction{T,P}},
) where {T,P}
    return MOI.VectorQuadraticFunction{T}
end
