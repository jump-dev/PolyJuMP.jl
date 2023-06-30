import DynamicPolynomials

const VariableOrder =
    DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}
const MonomialOrder = MP.Graded{MP.LexOrder}
const VarType = DynamicPolynomials.Variable{VariableOrder,MonomialOrder}
const PolyType{T} = DynamicPolynomials.Polynomial{VariableOrder,MonomialOrder,T}
const FuncType{T} = PolyJuMP.ScalarPolynomialFunction{T,PolyType{T}}

"""
    ToPolynomialBridge{T,S} <: Bridges.Constraint.AbstractBridge

`ToPolynomialBridge` implements the following reformulations:

  * ``f(x) \\in S`` into ``p(x) \\in S`` where ``f(x)`` is a scalar function
    and ``p(x)`` is a polynomial.

## Source node

`ToPolynomialBridge` supports:

  * `F` in `S` where `F` is a `MOI.AbstractScalarFunction` for which
   `convert(::Type{PolyJuMP.ScalarPolynomialFunction}, ::Type{F})`. That is for
   instance the case for `MOI.VariableIndex`, `MOI.ScalarAffineFunction` and
    `MOI.ScalarQuadraticFunction`.

## Target nodes

`ToPolynomialBridge` creates:

  * [`PolyJuMP.ScalarPolynomialFunction{T}`](@ref) in `S`
"""
struct ToPolynomialBridge{T,S} <:
       MOI.Bridges.Constraint.AbstractFunctionConversionBridge{FuncType{T},S}
    constraint::MOI.ConstraintIndex{FuncType{T},S}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ToPolynomialBridge{T,S}},
    model,
    f::MOI.AbstractScalarFunction,
    s::S,
) where {T,S}
    constraint = MOI.add_constraint(model, convert(FuncType{T}, f), s)
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

function MOI.modify(
    model::MOI.ModelLike,
    b::ToPolynomialBridge,
    change::MOI.AbstractFunctionModification,
)
    MOI.modify(model, b.constraint, change)
    return
end
