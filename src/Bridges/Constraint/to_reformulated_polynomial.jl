"""
    ToReformulatedPolynomialBridge{T,S} <: Bridges.Constraint.AbstractBridge

`ToReformulatedPolynomialBridge` implements a sparsity-preserving
reformulation of `f(x) in S` into `g(x, y) in S` where `f` is a
`MOI.ScalarNonlinearFunction`: products and powers of sums are replaced by
auxiliary variables `y`, so that converting `g` into a polynomial does not
expand them.

Since source and target are both `MOI.ScalarNonlinearFunction`, this bridge
must not be added to a `MOI.Bridges.LazyBridgeOptimizer` (it would bridge
its own output); use it in a `MOI.Bridges.Constraint.SingleBridgeOptimizer`
on top of an optimizer converting `MOI.ScalarNonlinearFunction` into
polynomials, e.g., [`PolyJuMP.QCQP.Optimizer`](@ref).

## Source node

`ToReformulatedPolynomialBridge` supports:

  * `MOI.ScalarNonlinearFunction` in `S`

## Target nodes

`ToReformulatedPolynomialBridge` creates:

  * `MOI.ScalarNonlinearFunction` in `S`
  * `MOI.ScalarNonlinearFunction` in `MOI.EqualTo{T}`, defining the
    auxiliary variables
  * `MOI.VariableIndex` in `MOI.Interval{T}`, the auxiliary variables
"""
struct ToReformulatedPolynomialBridge{T,S} <:
       MOI.Bridges.Constraint.AbstractFunctionConversionBridge{MOI.ScalarNonlinearFunction,S}
    constraint::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ToReformulatedPolynomialBridge{T,S}},
    model::MOI.ModelLike,
    func::MOI.ScalarNonlinearFunction,
    set::S,
) where {T,S}
    reformulated = PolyJuMP.Bridges._to_reformulated_polynomial!(model, T, func)
    constraint = MOI.add_constraint(model, reformulated, set)
    return ToReformulatedPolynomialBridge{T,S}(constraint)
end

function MOI.Bridges.Objective.supports_objective_function(
    ::Type{ToReformulatedPolynomialBridge{T}},
    ::Type{F},
) where {T,F<:MOI.ScalarNonlinearFunction}
    return MOI.Utilities.is_coefficient_type(F, T)
end


function MOI.supports_constraint(
    ::Type{ToReformulatedPolynomialBridge{T}},
    ::Type{<:MOI.ScalarNonlinearFunction},
    ::Type{<:MOI.AbstractScalarSet},
) where {T}
    return true
end


function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.Bridges.ObjectiveFunctionValue,
    ::ToReformulatedPolynomialBridge{T},
) where {T}
    F = FuncType{T}
    attr_f = MOI.Bridges.ObjectiveFunctionValue{F}(attr.result_index)
    return MOI.get(model, attr_f)
end
