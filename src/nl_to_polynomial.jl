# This will be refactored into a constraint bridge once https://github.com/jump-dev/MathOptInterface.jl/issues/846 is done

import DynamicPolynomials
const VarType = DynamicPolynomials.PolyVar{true}
const PolyType{T} = DynamicPolynomials.Polynomial{true,T}
const FuncType{T} = ScalarPolynomialFunction{T,PolyType{T}}

mutable struct NLToPolynomial{T,M<:MOI.ModelLike} <: MOI.AbstractOptimizer
    model::M
    constraint_indices::Vector{MOI.ConstraintIndex{FuncType{T}}}
    invalid::Union{Nothing,String}
    function NLToPolynomial{T}(model::M) where {T,M}
        return new{T,M}(model, MOI.ConstraintIndex{FuncType{T}}[], nothing)
    end
end

# The interesting part: rewriting NL into polynomial

struct InvalidNLExpression <: Exception
    message::String
end

function _to_polynomial!(d, ::Type, expr)
    throw(
        InvalidNLExpression(
            "Unexpected expression type `$(typeof(expr))`` of `$expr`",
        ),
    )
end

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
        throw(InvalidNLExpression("Unrecognized expression `$(expr)`"))
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

function _to_polynomial(model::NLToPolynomial{T}, expr) where {T}
    try
        return _to_polynomial(expr, T)
    catch err
        if err isa InvalidNLExpression
            model.invalid = string(
                "Cannot convert expression `$(expr)` into a polynomial: ",
                err.message,
                ".",
            )
            return
        else
            rethrow(err)
        end
    end
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
    model.invalid = nothing
    if data.has_objective
        obj = _to_polynomial(model, MOI.objective_expr(data.evaluator))
        if isnothing(obj)
            return
        end
        MOI.set(model.model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    end
    model.constraint_indices = map(eachindex(data.constraint_bounds)) do i
        func, set = MOI.FileFormats.MOF.extract_function_and_set(
            MOI.constraint_expr(data.evaluator, i),
        )
        poly = _to_polynomial(model, func)
        if isnothing(poly)
            return
        end
        return MOI.add_constraint(model, poly, set)
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
    MOI.empty!(model.model)
    model.invalid = nothing
    return
end

function MOI.is_empty(model::NLToPolynomial)
    return MOI.is_empty(model.model) && isnothing(model.invalid)
end

function MOI.optimize!(model::NLToPolynomial)
    if isnothing(model.invalid)
        MOI.optimize!(model.model)
    end
    return
end

_invalid_value(model::NLToPolynomial, ::MOI.RawStatusString) = model.invalid
_invalid_value(::NLToPolynomial, ::MOI.TerminationStatus) = MOI.INVALID_MODEL
_invalid_value(::NLToPolynomial, ::MOI.ResultCount) = 0
function _invalid_value(
    ::NLToPolynomial,
    ::Union{MOI.PrimalStatus,MOI.DualStatus},
)
    return MOI.NO_SOLUTION
end
function _invalid_value(
    ::NLToPolynomial,
    attr::Union{MOI.VariablePrimal,MOI.ConstraintDual,MOI.ConstraintPrimal},
)
    throw(MOI.ResultIndexBoundsError(attr, 0))
end

function MOI.get(model::NLToPolynomial, attr::MOI.AbstractModelAttribute)
    if isnothing(model.invalid) || !MOI.is_set_by_optimize(attr)
        return MOI.get(model.model, attr)
    else
        return _invalid_value(model, attr)
    end
end

function MOI.get(
    model::NLToPolynomial,
    attr::MOI.AbstractVariableAttribute,
    vi::MOI.VariableIndex,
)
    if isnothing(model.invalid) || !MOI.is_set_by_optimize(attr)
        return MOI.get(model.model, attr, vi)
    else
        return _invalid_value(model, attr)
    end
end

function MOI.get(
    model::NLToPolynomial,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
)
    if isnothing(model.invalid) || !MOI.is_set_by_optimize(attr)
        return MOI.get(model.model, attr, ci)
    else
        return _invalid_value(model, attr)
    end
end

# The boilerplate part: passing everything to the inner `.model`

function MOI.supports_incremental_interface(model::NLToPolynomial)
    return MOI.supports_incremental_interface(model.model)
end
function MOI.copy_to(dest::NLToPolynomial, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

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

function MOI.get(model::NLToPolynomial, attr::MOI.AbstractOptimizerAttribute)
    return MOI.get(model.model, attr)
end

function MOI.set(
    model::NLToPolynomial,
    attr::Union{MOI.AbstractModelAttribute,MOI.AbstractOptimizerAttribute},
    value,
)
    return MOI.set(model.model, attr, value)
end

function MOI.supports(
    model::NLToPolynomial,
    attr::MOI.AbstractVariableAttribute,
    ::Type{MOI.VariableIndex},
)
    return MOI.supports(model.model, attr, MOI.VariableIndex)
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

function MOI.set(
    model::NLToPolynomial,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
    value,
)
    return MOI.set(model.model, attr, ci, value)
end