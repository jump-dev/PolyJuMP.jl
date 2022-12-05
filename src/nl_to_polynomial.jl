# This will be moved to a bridge once https://github.com/jump-dev/MathOptInterface.jl/issues/846 is done

import DynamicPolynomials
const VarType = DynamicPolynomials.PolyVar{true}
const PolyType{T} = DynamicPolynomials.Polynomial{true,T}
const FuncType{T} = ScalarPolynomialFunction{T,PolyType{T}}

mutable struct NLToPolynomial{T,M<:MOI.ModelLike} <: MOI.AbstractOptimizer
    model::M
    constraint_indices::Vector{MOI.ConstraintIndex{FuncType{T}}}
    function NLToPolynomial{T}(model::M) where {T,M}
        return new{T,M}(model, MOI.ConstraintIndex{FuncType{T}}[])
    end
end

# The interesting part: rewriting NL into polynomial

_to_polynomial!(d, ::Type{T}, x::T) where {T} = x

function _to_polynomial!(d, ::Type{T}, expr::Expr) where {T}
    if Base.isexpr(expr, :call) && expr.args[1] === :+
        return sum(
            _to_polynomial!.(Ref(d), T, expr.args[2:end]),
            init = zero(PolyType{T}),
        )
    elseif Base.isexpr(expr, :call) && expr.args[1] === :*
        return prod(
            _to_polynomial!.(Ref(d), T, expr.args[2:end]),
            init = one(PolyType{T}),
        )
    elseif Base.isexpr(expr, :call) &&
           expr.args[1] === :- &&
           length(expr.args) == 2
        return -_to_polynomial!(d, T, expr.args[2])
    elseif Base.isexpr(expr, :call) &&
           expr.args[1] === :- &&
           length(expr.args) == 3
        a, b = _to_polynomial!.(Ref(d), T, expr.args[2:end])
        return a - b
    elseif Base.isexpr(expr, :ref) &&
           expr.args[1] === :x &&
           expr.args[2] isa MOI.VariableIndex
        vi = expr.args[2]
        if !haskey(d, vi)
            d[vi] = MP.similarvariable(VarType, Symbol("x[$(vi.value)]"))
        end
        return d[vi]
    else
        error("Cannot convert expression \"$(expr)\" into a polynomial.")
    end
end

function _to_polynomial(expr, ::Type{T}) where {T}
    d = Dict{MOI.VariableIndex,VarType}()
    poly = _to_polynomial!(d, T, expr)
    variable_map = collect(d)
    sort!(variable_map, by = x -> x[2])
    variables = [x[1] for x in variable_map]
    return FuncType{T}(poly, variables)
end

MOI.supports(::NLToPolynomial, ::MOI.NLPBlock) = true
function MOI.set(
    model::NLToPolynomial{T},
    ::MOI.NLPBlock,
    data::MOI.NLPBlockData,
) where {T}
    # FIXME if a non-NLP objective is set afterwards, it might overwrite.
    # but let's not complicate as it will be fixed by
    # https://github.com/jump-dev/MathOptInterface.jl/issues/846
    MOI.initialize(data.evaluator, [:ExprGraph])
    if data.has_objective
        obj = _to_polynomial(MOI.objective_expr(data.evaluator), T)
        MOI.set(model.model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    end
    model.constraint_indices = map(eachindex(data.constraint_bounds)) do i
        func, set = MOI.FileFormats.MOF.extract_function_and_set(
            MOI.constraint_expr(data.evaluator, i),
        )
        return MOI.add_constraint(model, _to_polynomial(func, T), set)
    end
end

function MOI.get(model::NLToPolynomial, attr::MOI.NLPBlockDual)
    return map(model.constraint_indices) do ci
        return MOI.get(model.model, MOI.ConstraintDual(attr.result_index), ci)
    end
end

# TODO not used yet and we don't want to use MOI.Nonlinear as it might break with minor releases
#function MOI.get(model::NLToPolynomial, attr::MOI.ConstraintDual, ci::MOI.Nonlinear.ConstraintIndex)
#    return MOI.get(model.model, attr, model.constraint_indices[ci.value])
#end

function MOI.empty!(model::NLToPolynomial)
    empty!(model.constraint_indices)
    return MOI.empty!(model.model)
end

# The boilerplate part: passing everything to the inner `.model`

MOI.is_empty(model::NLToPolynomial) = MOI.is_empty(model.model)

function MOI.supports_incremental_interface(model::NLToPolynomial)
    return MOI.supports_incremental_interface(model.model)
end
function MOI.copy_to(dest::NLToPolynomial, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

MOI.optimize!(model::NLToPolynomial) = MOI.optimize!(model.model)

MOI.add_variable(model::NLToPolynomial) = MOI.add_variable(model.model)

function MOI.supports_constraint(
    model::NLToPolynomial,
    ::Type{F},
    ::Type{S},
) where {F<:MOI.AbstractFunction,S<:MOI.AbstractSet}
    return MOI.supports_constraint(model.model, F, S)
end

function MOI.add_constraint(
    model::NLToPolynomial,
    func::MOI.AbstractFunction,
    set::MOI.AbstractSet,
)
    return MOI.add_constraint(model.model, func, set)
end

function MOI.supports(
    model::NLToPolynomial,
    attr::Union{MOI.AbstractOptimizerAttribute,MOI.AbstractModelAttribute},
)
    return MOI.supports(model.model, attr)
end

function MOI.get(
    model::NLToPolynomial,
    attr::Union{MOI.AbstractOptimizerAttribute,MOI.AbstractModelAttribute},
)
    return MOI.get(model.model, attr)
end

function MOI.set(model::NLToPolynomial, attr::MOI.AbstractModelAttribute, value)
    return MOI.set(model.model, attr, value)
end

function MOI.supports(
    model::NLToPolynomial,
    attr::MOI.AbstractVariableAttribute,
    ::Type{MOI.VariableIndex},
)
    return MOI.supports(model.model, attr, MOI.VariableIndex)
end

function MOI.get(
    model::NLToPolynomial,
    attr::MOI.AbstractVariableAttribute,
    vi::MOI.VariableIndex,
)
    return MOI.get(model.model, attr, vi)
end

function MOI.set(
    model::NLToPolynomial,
    attr::MOI.AbstractVariableAttribute,
    vi::MOI.VariableIndex,
    value,
)
    return MOI.set(model.model, attr, vi, value)
end

function MOI.supports(
    model::NLToPolynomial,
    attr::MOI.AbstractConstraintAttribute,
    ::Type{MOI.ConstraintIndex{F,S}},
) where {F,S}
    return MOI.supports(model.model, attr, MOI.ConstraintIndex{F,S})
end

function MOI.get(
    model::NLToPolynomial,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
)
    return MOI.get(model.model, attr, ci)
end

function MOI.set(
    model::NLToPolynomial,
    attr::MOI.AbstractConstraintAttribute,
    vi::MOI.ConstraintIndex,
    value,
)
    return MOI.set(model.model, attr, ci, value)
end
