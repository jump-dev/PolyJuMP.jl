import MultivariatePolynomials as MP
import DynamicPolynomials


struct ToReformulatedPolynomialBridge{T} <: MOI.Bridges.Objective.AbstractBridge end

function MOI.Bridges.Objective.bridge_objective(
    ::Type{ToReformulatedPolynomialBridge{T}},
    model::MOI.ModelLike,
    func::MOI.ScalarNonlinearFunction,
) where {T}
    F = FuncType{T}
	nonlin_func = convert(MOI.ScalarNonlinearFunction, func)
	poly = PolyJuMP.Bridges._to_reformulated_polynomial!(model, T, nonlin_func)
	MOI.set(model, MOI.ObjectiveFunction{F}(), poly)
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
