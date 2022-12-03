import DynamicPolynomials

const VarType = DynamicPolynomials.PolyVar{true}
const PolyType{T} = DynamicPolynomials.Polynomial{true,T}
const FuncType{T} = ScalarPolynomialFunction{T,DynamicPolynomials.Polynomial{true,T}}

"""
    ToPolynomialBridge{T,S} <: Bridges.Constraint.AbstractBridge

`ToPolynomialBridge` implements the following reformulations:

  * ``f(x) \\in S`` into ``p(x) \\in S`` where ``f(x)`` is either
    affine or quadratic and ``p(x)`` is a polynomial.`

## Source node

`ToPolynomialBridge` supports:

  * [`MOI.VariableIndex`](@ref) in `S`
  * [`MOI.ScalarAffineFunction`](@ref) in `S`
  * [`MOI.ScalarQuadraticFunction`](@ref) in `S`

## Target nodes

`ToPolynomialBridge` creates:

  * [`MOI.ScalarPolynomialFunction{T}`](@ref) in `S`
"""
struct ToPolynomialBridge{T,S} <:
       MOI.Bridges.Constraint.AbstractFunctionConversionBridge{FuncType{T},S}
    constraint::MOI.ConstraintIndex{FuncType{T},S}
end

function _to_polynomial(::Type{T}, vi::MOI.VariableIndex) where T
    DynamicPolynomials.@polyvar x
    return FuncType{T}(polynomial(x, T), [vi])
end

function _to_polynomial(::Type{T}, func::MOI.ScalarAffineFunction{T}) where T
    variables = [t.variable for t in func.terms]
    sort!(variables)
    unique!(variables)
    DynamicPolynomials.@polyvar x[1:length(variables)]
    d = Dict(variables[i] => x[i] for i in eachindex(variables))
    terms = MP.termtype(VarType, T)[MOI.constant(func)]
    for t in func.terms
        push!(terms, MP.term(t.coefficient, d[t.variable]))
    end
    return FuncType{T}(polynomial(terms), variables)
end

function _to_polynomial(::Type{T}, func::MOI.ScalarQuadraticFunction{T}) where T
    linear_variables = [t.variable for t in func.affine_terms]
    quad_variables_1 = [t.variable_1 for t in func.quadratic_terms]
    quad_variables_2 = [t.variable_2 for t in func.quadratic_terms]
    variables = [linear_variables; quad_variables_1; quad_variables_2]
    sort!(variables, by = v -> v.value)
    unique!(variables)
    DynamicPolynomials.@polyvar x[1:length(variables)]
    d = Dict(variables[i] => x[i] for i in eachindex(variables))
    terms = MP.termtype(VarType, T)[MOI.constant(func)]
    for t in func.affine_terms
        push!(terms, MP.term(t.coefficient, d[t.variable]))
    end
    for t in func.quadratic_terms
        coef = t.variable_1 == t.variable_2 ? t.coefficient / 2 : t.coefficient
        push!(terms, MP.term(coef, d[t.variable_1] * d[t.variable_2]))
    end
    return FuncType{T}(polynomial(terms), variables)
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ToPolynomialBridge{T,S}},
    model,
    f::MOI.AbstractScalarFunction,
    s::S,
) where {T,S}
    constraint = MOI.add_constraint(model, _to_polynomial(T, f), s)
    return ToPolynomialBridge{T,S}(constraint)
end

function MOI.supports_constraint(
    ::Type{ToPolynomialBridge{T}},
    ::Type{<:MOI.AbstractScalarFunction},
    ::Type{<:MOI.AbstractScalarSet},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:ToPolynomialBridge},
)
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{ToPolynomialBridge{T,S}},
) where {T,S}
    return Tuple{Type,Type}[(FuncType{T}, S)]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ToPolynomialBridge{T}},
    ::Type{<:MOI.AbstractScalarFunction},
    ::Type{S},
) where {T,S<:MOI.AbstractScalarSet}
    return ToPolynomialBridge{T,S}
end

function MOI.get(
    ::ToPolynomialBridge{T,S},
    ::MOI.NumberOfConstraints{FuncType{T},S},
)::Int64 where {T,S}
    return 1
end

function MOI.get(
    b::ToPolynomialBridge{T,S},
    ::MOI.ListOfConstraintIndices{FuncType{T},S},
) where {T,S}
    return [b.constraint]
end

function MOI.delete(model::MOI.ModelLike, b::ToPolynomialBridge)
    MOI.delete(model, b.constraint)
    return
end
