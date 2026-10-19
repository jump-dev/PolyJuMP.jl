import MultivariatePolynomials as MP
import DynamicPolynomials
using IntervalArithmetic
const IA = IntervalArithmetic

const VariableOrder =
    DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}
const MonomialOrder = MP.Graded{MP.LexOrder}
const VarType = DynamicPolynomials.Variable{VariableOrder,MonomialOrder}
const PolyType{T} = DynamicPolynomials.Polynomial{VariableOrder,MonomialOrder,T}
const FuncType{T} = PolyJuMP.ScalarPolynomialFunction{T,PolyType{T}}

# Auxiliary variables created by the sparsity-preserving reformulation, keyed
# by `_expr_key` of the expression they replace.
const ExprMap = Dict{Any,MOI.VariableIndex}

_expr_key(x::MOI.VariableIndex) = x
_expr_key(x::Number) = x
_expr_key(x) = _expr_key(convert(MOI.ScalarNonlinearFunction, x))
function _expr_key(f::MOI.ScalarNonlinearFunction)
    args = Any[_expr_key(arg) for arg in f.args]
    if f.head in (:+, :*)
        sort!(args; by = string)
    end
    return (f.head, args...)
end


"""
    _to_reformulated_polynomial!(model, T, expr)

Returns `expr` where products and powers of sums whose expansion would be
large are replaced by auxiliary variables :

* in a product of degree larger than 2 with several factors containing a
  sum, these factors are substituted from the highest degree to the lowest
  until only one remains, since expanding the product of one sum with
  monomials does not increase the density much;
* the base of a power is substituted unless the power has degree less than
  3, or degree 3 with at most 2 variables.

The auxiliary variables are added to `model` with bounds computed by
interval arithmetic and their defining equality constraints.
"""
function _to_reformulated_polynomial!(
    model::MOI.ModelLike,
    ::Type{T},
    expr,
) where {T}
    d = ExprMap()
    return _to_reformulated_polynomial!(model, T, d, expr)
end

function _to_reformulated_polynomial!(
    model::MOI.ModelLike,
    ::Type{T},
    d::ExprMap,
    expr,
) where {T}
    return _to_reformulated_polynomial!(model, T, d, convert(MOI.ScalarNonlinearFunction, expr))
end

function _to_reformulated_polynomial!(
    ::MOI.ModelLike,
    ::Type{T},
    ::ExprMap,
    expr::Union{MOI.VariableIndex,Number},
) where {T}
    return expr
end

function _to_reformulated_polynomial!(
    model::MOI.ModelLike,
    ::Type{T},
    d::ExprMap,
    expr::MOI.ScalarNonlinearFunction,
) where {T}
    operator = expr.head
    operands = expr.args
    if operator == :*
        remaining = count(_is_expr, operands)
        if _degree(expr) > 2 && remaining > 1
            new_operands = Vector{Any}(undef, length(operands))
            for i in sortperm(operands; by = _degree, rev = true)
                op = operands[i]
                if remaining > 1 && _is_expr(op)
                    new_operands[i] = _substitute_expr(model, T, d, op)
                    remaining -= 1
                else
                    new_operands[i] = _to_reformulated_polynomial!(model, T, d, op)
                end
            end
            return MOI.ScalarNonlinearFunction(:*, new_operands)
        end
    elseif operator == :^ && length(operands) == 2
        base, exponent = operands
        deg = _degree(expr)
        if deg < 3 || (deg == 3 && length(_var_index_set(expr)) <= 2)
            base = _to_reformulated_polynomial!(model, T, d, base)
        else
            base = _substitute_expr(model, T, d, base)
        end
        return MOI.ScalarNonlinearFunction(:^, Any[base, exponent])
    end
    return MOI.ScalarNonlinearFunction(
        operator,
        Any[_to_reformulated_polynomial!(model, T, d, op) for op in operands],
    )
end

_degree(::Number) = 0
_degree(::MOI.VariableIndex) = 1
_degree(expr) = _degree(convert(MOI.ScalarNonlinearFunction, expr))

function _degree(expr::MOI.ScalarNonlinearFunction)
    operator = expr.head
    operands = expr.args

    if operator == :+ || operator == :-
        return maximum(_degree, operands)
    elseif operator == :*
        return sum(_degree, operands)
    elseif operator == :/ && length(operands) == 2
        return _degree(operands[1])
    elseif operator == :^ && length(operands) == 2
        if _degree(operands[2]) != 0
            error("Cannot infer the degree of $expr: the exponent is not constant.")
        end
        return _degree(operands[1]) * Int(_evaluate_const_expr(operands[2]))
    end
    return error("Unknown expression $expr, cannot infer degree")
end

_evaluate_const_expr(expr::MOI.VariableIndex) = error("Cannot evaluate $expr, expr needs to be constant")
_evaluate_const_expr(expr::Number) = expr
_evaluate_const_expr(expr) = _evaluate_const_expr(convert(MOI.ScalarNonlinearFunction, expr))

function _evaluate_const_expr(expr::MOI.ScalarNonlinearFunction)
    operator = expr.head
    operands = expr.args
    if operator == :+
        return sum(_evaluate_const_expr, operands)
    elseif operator == :- && length(operands) == 2
        return _evaluate_const_expr(operands[1]) - _evaluate_const_expr(operands[2])
    elseif operator == :- && length(operands) == 1
        return -_evaluate_const_expr(operands[1])
    elseif operator == :*
        return prod(_evaluate_const_expr, operands)
    elseif operator == :/ && length(operands) == 2
        return _evaluate_const_expr(operands[1]) / _evaluate_const_expr(operands[2])
    elseif operator == :^ && length(operands) == 2
        return _evaluate_const_expr(operands[1])^_evaluate_const_expr(operands[2])
    end
    return error("Cannot evaluate $expr")
end

_var_index_set(::Number, ::Set{Int}) = nothing
_var_index_set(expr::MOI.VariableIndex, vars::Set{Int}) = push!(vars, expr.value)
function _var_index_set(expr, vars::Set{Int})
    return _var_index_set(convert(MOI.ScalarNonlinearFunction, expr), vars)
end

function _var_index_set(expr)::Set{Int}
    ret = Set{Int}()
    return _var_index_set(convert(MOI.ScalarNonlinearFunction, expr), ret)
end

function _var_index_set(expr::MOI.ScalarNonlinearFunction, vars::Set{Int})
    for operand in expr.args
        _var_index_set(operand, vars)
    end
    return vars
end

"""
    _substitute_expr(model, T, d::ExprMap, expr)

If `expr` is a sum or a difference, returns an auxiliary variable equal to
its reformulation (reusing the one created for an identical expression of
the same function); otherwise returns the reformulation of `expr`.
"""
function _substitute_expr(
    model::MOI.ModelLike,
    ::Type{T},
    d::ExprMap,
    expr,
) where {T}
    func = convert(MOI.ScalarNonlinearFunction, expr)
    func = _to_reformulated_polynomial!(model, T, d, func)
    if !(func.head in (:+, :-))
        return func
    end
    key = _expr_key(func)
    v = get(d, key, nothing)
    if v === nothing
        lb, ub = IA.bounds(_get_nonlinfunc_bounds(model, T, func))
        v, _ = MOI.add_constrained_variable(model, MOI.Interval(lb, ub))
        MOI.Utilities.normalize_and_add_constraint(
            model,
            MOI.ScalarNonlinearFunction(:-, Any[func, v]),
            MOI.EqualTo(zero(T)),
        )
        d[key] = v
    end
    return v
end

function _substitute_expr(
    ::MOI.ModelLike,
    ::Type{T},
    ::ExprMap,
    expr::Union{MOI.VariableIndex,Number},
) where {T}
    return expr
end

"""
    _is_expr(expr)

Returns whether `expr` contains a sum or a difference of several terms.
"""
function _is_expr(expr::MOI.ScalarNonlinearFunction)::Bool
    return (expr.head in (:+, :-) && length(expr.args) > 1) || any(_is_expr, expr.args)
end
_is_expr(::Union{MOI.VariableIndex,Number}) = false
_is_expr(expr) = _is_expr(convert(MOI.ScalarNonlinearFunction, expr))

"""
    _get_nonlinfunc_bounds(model, T, expr)

Interval containing the values of `expr` given the bounds of the variables
in `model`, used as bounds of the auxiliary variable replacing `expr`. A
wrong (too narrow) interval cuts feasible solutions, so every operation must
return a valid enclosure.
"""
function _get_nonlinfunc_bounds(model::MOI.ModelLike, ::Type{T}, expr) where {T}
    return _get_nonlinfunc_bounds(
        model,
        T,
        convert(MOI.ScalarNonlinearFunction, expr),
    )
end

function _get_nonlinfunc_bounds(::MOI.ModelLike, ::Type{T}, expr::Number) where {T}
    return IA.interval(float(expr))
end

function _get_nonlinfunc_bounds(
    model::MOI.ModelLike,
    ::Type{T},
    expr::MOI.VariableIndex,
) where {T}
    F = float(T)
    lb, ub = MOI.Utilities.get_bounds(model, T, expr)
    return IA.interval(
        lb == typemin(lb) ? typemin(F) : float(lb),
        ub == typemax(ub) ? typemax(F) : float(ub),
    )
end

function _get_nonlinfunc_bounds(
    model::MOI.ModelLike,
    ::Type{T},
    expr::MOI.ScalarNonlinearFunction,
) where {T}
    operator = expr.head
    operands = expr.args
    bounds(op) = _get_nonlinfunc_bounds(model, T, op)
    if operator == :+
        return sum(bounds, operands)
    elseif operator == :- && length(operands) == 1
        return -bounds(operands[1])
    elseif operator == :- && length(operands) == 2
        return bounds(operands[1]) - bounds(operands[2])
    elseif operator == :*
        return prod(bounds, operands)
    elseif operator == :/ && length(operands) == 2
        return bounds(operands[1]) / bounds(operands[2])
    elseif operator == :^ && length(operands) == 2
        base, exponent = operands
        if _degree(exponent) == 0 && isinteger(_evaluate_const_expr(exponent))
            return bounds(base)^Int(_evaluate_const_expr(exponent))
        end
        return bounds(base)^bounds(exponent)
    end
    return error("Unknown expression $expr, cannot compute its bounds")
end
