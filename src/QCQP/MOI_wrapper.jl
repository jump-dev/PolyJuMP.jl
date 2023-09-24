import MathOptInterface as MOI

mutable struct Optimizer{T,O<:MOI.ModelLike} <: MOI.AbstractOptimizer
    model::O
    objective::Union{Nothing,PolyJuMP.ScalarPolynomialFunction{T}}
end

function Optimizer(model::MOI.ModelLike)
    return Optimizer{Float64,typeof(model)}(model, nothing)
end

function MOI.get(
    model::Optimizer{T},
    attr::MOI.Bridges.ListOfNonstandardBridges,
) where {T}
    list = copy(MOI.get(model.model, attr))
    push!(list, PolyJuMP.Bridges.Constraint.ToPolynomialBridge{T})
    push!(list, PolyJuMP.Bridges.Objective.ToPolynomialBridge{T})
    return list
end

MOI.is_empty(model::Optimizer) = MOI.is_empty(model.model)
MOI.empty!(model::Optimizer) = MOI.empty!(model.model)

MOI.add_variable(model::Optimizer) = MOI.add_variable(model.model)

function MOI.supports_add_constrained_variable(
    model::Optimizer,
    ::Type{S},
) where {S<:MOI.AbstractScalarSet}
    return MOI.supports_add_constrained_variable(model.model, S)
end

function MOI.supports_add_constrained_variables(
    model::Optimizer,
    ::Type{MOI.Reals},
)
    return MOI.supports_add_constrained_variables(model.model, MOI.Reals)
end

function MOI.supports_add_constrained_variables(
    model::Optimizer,
    ::Type{S},
) where {S<:MOI.AbstractVectorSet}
    return MOI.supports_add_constrained_variables(model.model, S)
end

function MOI.supports(model::Optimizer, attr::MOI.AbstractModelAttribute)
    return MOI.supports(model.model, attr)
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{<:PolyJuMP.ScalarPolynomialFunction},
)
    return true
end

function MOI.set(
    model::Optimizer{T},
    ::MOI.ObjectiveFunction{F},
    f::F,
) where {T,F<:PolyJuMP.ScalarPolynomialFunction{T}}
    model.objective = f
    return
end

function MOI.set(model::Optimizer, attr::MOI.AbstractModelAttribute, value)
    MOI.set(model.model, attr, value)
end

function MOI.get(model::Optimizer, attr::MOI.AbstractModelAttribute)
    return MOI.get(model.model, attr)
end


function MOI.supports_constraint(
    model::Optimizer,
    ::Type{F},
    ::Type{S},
) where {F<:MOI.AbstractFunction,S<:MOI.AbstractSet}
    return MOI.supports_constraint(model.model, F, S)
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.AbstractFunction,
    set::MOI.AbstractSet,
)
    return MOI.add_constraint(model.model, func, set)
end

function MOI.supports_incremental_interface(model::Optimizer)
    return MOI.supports_incremental_interface(model.model)
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

MOI.optimize!(model::Optimizer) = MOI.optimize!(model.model)

function _quad_convert(index, p::MP.AbstractPolynomialLike{T}, div) where {T}
    q = zero(MOI.ScalarQuadraticFunction{T})
    for t in MP.terms(p)
        α = MP.coefficient(t)
        mono = MP.monomial(t)
        if MP.degree(mono) == 0
            MA.operate!(+, q, α)
        else
            vars = MP.effective_variables(mono)
            if MP.degree(mono) == 1
                @assert length(vars) == 1
                MA.operate!(MA.add_mul, q, α, index(first(vars)))
            elseif MP.degree(mono) == 2
                x = first(vars)
                if length(vars) == 1
                    y = x
                else
                    @assert length(vars) == 2
                    y = vars[2]
                end
                MA.operate!(MA.add_mul, q, α, index(x), index(y))
            else
                x = div[mono]
                y = MP.div_multiple(mono, x)
                MA.operate!(MA.add_mul, q, α, index(x), index(y))
            end
        end
    end
    return q
end

function _variable_to_index_map(p::PolyJuMP.ScalarPolynomialFunction{T,P}) where {T,P}
    M = MP.monomial_type(P)
    return Dict{M,MOI.VariableIndex}(
        v => vi for (v, vi) in zip(MP.variables(p.polynomial), p.variables)
    )
end

function monomial_variable_index(model::Optimizer{T}, d::Dict, div, mono::MP.AbstractMonomialLike) where {T}
    if !haskey(d, mono)
        d[mono] = MOI.add_variable(model.model)
        x = div[mono]
        vx = monomial_variable_index(model, d, div, x)
        y = MP.div_multiple(mono, x)
        vy = monomial_variable_index(model, d, div, y)
        MOI.add_constraint(model,
            MA.@rewrite(one(T) * d[mono] - one(T) * vx * vy),
            MOI.EqualTo(zero(T)),
        )
    end
    return d[mono]
end

function MOI.Utilities.final_touch(model::Optimizer{T}, _) where {T}
    if !isnothing(model.objective)
        d = _variable_to_index_map(model.objective)
        p = model.objective.polynomial
        monos = MP.monomials(p)
        div = decompose(monos)
        function index(mono)
            return monomial_variable_index(model, d, div, mono)
        end
        obj = _quad_convert(index, p, div)
        MOI.set(model.model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    end
    return
end

function MOI.get(model::Optimizer, attr::MOI.SolverName)
    name = MOI.get(model.model, attr)
    return "PolyJuMP.QCQP with $name"
end
