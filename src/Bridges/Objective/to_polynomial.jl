import DynamicPolynomials

const VarType = DynamicPolynomials.PolyVar{true}
const PolyType{T} = DynamicPolynomials.Polynomial{true,T}
const FuncType{T} = PolyJuMP.ScalarPolynomialFunction{T,PolyType{T}}

"""
    ToPolynomialBridge{T}

`ToPolynomialBridge` implements the following reformulations:

 * ``\\min \\{f(x)\\}`` into ``\\min\\{p(x)\\}``
 * ``\\max \\{f(x)\\}`` into ``\\max\\{p(x)\\}``

where ``f(x)`` is a scalar function and ``p(x)`` is a polynomial.

## Source node

`ToPolynomialBridge` supports:

 * `MOI.ObjectiveFunction{F}` where `F` is a `MOI.AbstractScalarFunction` for
   which `convert(::Type{PolyJuMP.ScalarPolynomialFunction}, ::Type{F})`. That
   is for instance the case for `MOI.VariableIndex`, `MOI.ScalarAffineFunction`
   and `MOI.ScalarQuadraticFunction`.

## Target nodes

`ToPolynomialBridge` creates:

 * One objective node: `MOI.ObjectiveFunction{PolyJuMP.ScalarPolynomialFunction{T}}`
"""
struct ToPolynomialBridge{T} <: MOI.Bridges.Objective.AbstractBridge end

function MOI.Bridges.Objective.bridge_objective(
    ::Type{ToPolynomialBridge{T}},
    model::MOI.ModelLike,
    func::MOI.AbstractScalarFunction,
) where {T}
    F = FuncType{T}
    MOI.set(model, MOI.ObjectiveFunction{F}(), convert(F, func))
    return ToPolynomialBridge{T}()
end

function MOI.Bridges.Objective.supports_objective_function(
    ::Type{ToPolynomialBridge{T}},
    ::Type{F},
) where {T,F<:MOI.AbstractScalarFunction}
    return MOI.Utilities.is_coefficient_type(F, T)
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:ToPolynomialBridge},
)
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(::Type{<:ToPolynomialBridge})
    return Tuple{Type,Type}[]
end

function MOI.Bridges.set_objective_function_type(
    ::Type{ToPolynomialBridge{T}},
) where {T}
    return FuncType{T}
end

# Attributes, Bridge acting as a model
MOI.get(::ToPolynomialBridge, ::MOI.NumberOfVariables)::Int64 = 0

function MOI.get(::ToPolynomialBridge, ::MOI.ListOfVariableIndices)
    return MOI.VariableIndex[]
end

# No variables or constraints are created in this bridge so there is nothing to
# delete.
MOI.delete(::MOI.ModelLike, ::ToPolynomialBridge) = nothing

function MOI.set(
    ::MOI.ModelLike,
    ::MOI.ObjectiveSense,
    ::ToPolynomialBridge,
    ::MOI.OptimizationSense,
)
    # `ToPolynomialBridge` is sense agnostic, therefore, we don't need to change
    # anything.
    return
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.Bridges.ObjectiveFunctionValue,
    ::ToPolynomialBridge{T},
) where {T}
    F = FuncType{T}
    attr_f = MOI.Bridges.ObjectiveFunctionValue{F}(attr.result_index)
    return MOI.get(model, attr_f)
end
