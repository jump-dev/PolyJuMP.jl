# This will be refactored into a constraint bridge once https://github.com/jump-dev/MathOptInterface.jl/issues/846 is done

import DynamicPolynomials
const VariableOrder =
    DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}
const MonomialOrder = MP.Graded{MP.LexOrder}
const VarType = DynamicPolynomials.Variable{VariableOrder,MonomialOrder}
const PolyType{T} = DynamicPolynomials.Polynomial{VariableOrder,MonomialOrder,T}
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
            "Unexpected expression type `$(typeof(expr))` of `$expr`",
        ),
    )
end

_to_polynomial!(d, ::Type{T}, x::T) where {T} = x
_to_polynomial!(d, ::Type, x::Number) = x

function _is_operator(expr::Expr, sym)
    return Base.Meta.isexpr(expr, :call) && expr.args[1] === sym
end
_operands(expr::Expr) = expr.args[2:end]
function _is_variable(expr::Expr)
    return Base.Meta.isexpr(expr, :ref) &&
           expr.args[1] === :x &&
           expr.args[2] isa MOI.VariableIndex
end

function _is_operator(func::MOI.ScalarNonlinearFunction, sym)
    return func.head === sym
end
_operands(func::MOI.ScalarNonlinearFunction) = func.args
_is_variable(expr::MOI.ScalarNonlinearFunction) = false

function _to_polynomial!(d, ::Type{T}, vi::MOI.VariableIndex) where {T}
    if !haskey(d, vi)
        d[vi] = MP.similar_variable(VarType, Symbol("x[$(vi.value)]"))
    end
    return d[vi]
end

function _to_polynomial!(
    d,
    ::Type{T},
    f::Union{MOI.ScalarAffineFunction,MOI.ScalarQuadraticFunction},
) where {T}
    return _to_polynomial!(d, T, convert(MOI.ScalarNonlinearFunction, f))
end

function _to_polynomial!(
    d,
    ::Type{T},
    expr::Union{Expr,MOI.ScalarNonlinearFunction},
) where {T}
    operands = _operands(expr)
    if _is_operator(expr, :+)
        return sum(
            _to_polynomial!.(Ref(d), T, operands),
            init = zero(PolyType{T}),
        )
    elseif _is_operator(expr, :*)
        return prod(
            _to_polynomial!.(Ref(d), T, operands),
            init = one(PolyType{T}),
        )
    elseif _is_operator(expr, :-) && length(operands) == 1
        return -_to_polynomial!(d, T, operands[1])
    elseif _is_operator(expr, :-) && length(operands) == 2
        a, b = _to_polynomial!.(Ref(d), T, operands)
        return a - b
    elseif _is_operator(expr, :^) && length(operands) == 2
        a, b = _to_polynomial!.(Ref(d), T, operands)
        if !(b isa Integer && round(Int, b) == b)
            b = round(Int, b)
        end
        return a^b
    elseif _is_operator(expr, :/) && length(operands) == 2
        a, b = _to_polynomial!.(Ref(d), T, operands)
        divisor, remainder = Base.divrem(a, b)
        if iszero(remainder)
            return divisor
        else
            return a / b
        end
    elseif _is_variable(expr)
        return _to_polynomial!(d, T, operands[1])
    else
        throw(InvalidNLExpression("Cannot convert `$(expr)` into a polynomial"))
    end
end

function _to_polynomial(expr, ::Type{T}) where {T}
    d = Dict{MOI.VariableIndex,VarType}()
    poly = _to_polynomial!(d, T, expr)
    return _scalar_polynomial(d, T, poly)
end

function _scalar_polynomial(d::Dict{K,V}, ::Type{T}, poly) where {T,K,V}
    inv = Dict(v => k for (k, v) in d)
    variables = [inv[v] for v in MP.variables(poly)]
    P = MP.polynomial_type(V, T)
    return ScalarPolynomialFunction{T,P}(poly, variables)
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
