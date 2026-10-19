"""
    ToReformulatedPolynomialBridge{T} <: Bridges.Objective.AbstractBridge

`ToReformulatedPolynomialBridge` implements a sparsity-preserving
reformulation of the objective `f(x)` into `g(x, y)` where `f` is a
`MOI.ScalarNonlinearFunction`: products and powers of sums are replaced by
auxiliary variables `y`, so that converting `g` into a polynomial does not
expand them.

Since source and target are both `MOI.ScalarNonlinearFunction`, this bridge
must not be added to a `MOI.Bridges.LazyBridgeOptimizer` (it would bridge
its own output); use it in a `MOI.Bridges.Objective.SingleBridgeOptimizer`
on top of an optimizer converting `MOI.ScalarNonlinearFunction` into
polynomials, e.g., [`PolyJuMP.QCQP.Optimizer`](@ref).

## Source node

`ToReformulatedPolynomialBridge` supports:

 * `MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}`

## Target nodes

`ToReformulatedPolynomialBridge` creates:

 * `MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}`
 * `MOI.ScalarNonlinearFunction` in `MOI.EqualTo{T}`, defining the auxiliary
   variables
 * `MOI.VariableIndex` in `MOI.Interval{T}`, the auxiliary variables
"""
struct ToReformulatedPolynomialBridge{T} <: MOI.Bridges.Objective.AbstractBridge end

function MOI.Bridges.Objective.bridge_objective(
    ::Type{ToReformulatedPolynomialBridge{T}},
    model::MOI.ModelLike,
    func::MOI.ScalarNonlinearFunction,
) where {T}
    reformulated = PolyJuMP.Bridges._to_reformulated_polynomial!(model, T, func)
    F = MOI.ScalarNonlinearFunction
    MOI.set(model, MOI.ObjectiveFunction{F}(), reformulated)
    return ToReformulatedPolynomialBridge{T}()
end

function MOI.Bridges.Objective.supports_objective_function(
    ::Type{ToReformulatedPolynomialBridge{T}},
    ::Type{F},
) where {T,F<:MOI.ScalarNonlinearFunction}
    return MOI.Utilities.is_coefficient_type(F, T)
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
