import MultivariatePolynomials as MP
import DynamicPolynomials


# struct ToReformulatedPolynomialBridge{T,S} <:
#        MOI.Bridges.Constraint.AbstractFunctionConversionBridge{FuncType{T},S}
#     constraint::MOI.ConstraintIndex{FuncType{T},S}
# end
struct ToReformulatedPolynomialBridge{T,S} <:
       MOI.Bridges.Constraint.AbstractFunctionConversionBridge{MOI.ScalarNonlinearFunction,S}
    constraint::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ToReformulatedPolynomialBridge{T}},
    model::MOI.ModelLike,
    func::MOI.ScalarNonlinearFunction,
	s::S,
) where {T,S}
	nonlin_func = convert(MOI.ScalarNonlinearFunction, func)
	poly = PolyJuMP.Bridges._to_reformulated_polynomial!(model, T, nonlin_func)
	# poly_func = convert(F, poly)
	# constraint = MOI.add_constraint(model, poly_func, s)
	constraint = MOI.add_constraint(model, poly, s)
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
