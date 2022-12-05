# This will be moved to a bridge once https://github.com/jump-dev/MathOptInterface.jl/issues/846 is done

import DynamicPolynomials
const VarType = DynamicPolynomials.PolyVar{true}
const PolyType{T} = DynamicPolynomials.Polynomial{true,T}
const FuncType{T} = ScalarPolynomialFunction{T,PolyType{T}}

struct NLToPolynomial{T,M<:MOI.ModelLike} <: MOI.AbstractOptimizer
    model::M
    function NLToPolynomial{T}(model::M) where {T,M}
        return new{T,M}(model)
    end
end

# The interesting part: rewriting NL into polynomial

_to_polynomial!(d, ::Type{T}, x::T) where {T} = x

function _to_polynomial!(d, ::Type{T}, expr::Expr) where {T}
    if Base.isexpr(expr, :call) && expr.args[1] === :+
        return sum(_to_polynomial!.(Ref(d), T, expr.args[2:end]), init=zero(PolyType{T}))
    elseif Base.isexpr(expr, :call) && expr.args[1] === :*
        return prod(_to_polynomial!.(Ref(d), T, expr.args[2:end]), init=one(PolyType{T}))
    elseif Base.isexpr(expr, :call) && expr.args[1] === :- && length(expr.args) == 2
        return -_to_polynomial!(d, T, expr.args[2])
    elseif Base.isexpr(expr, :call) && expr.args[1] === :- && length(expr.args) == 3
        a, b = _to_polynomial!.(Ref(d), T, expr.args[2:end])
        return a - b
    elseif Base.isexpr(expr, :ref) && expr.args[1] === :x && expr.args[2] isa MOI.VariableIndex
        vi = expr.args[2]
        if !haskey(d, vi)
            d[vi] = MP.similarvariable(VarType, Symbol("x[$(vi.value)]"))
        end
        return d[vi]
    else
        error("Cannot convert expression \"$(expr)\" into a polynomial.")
    end
end

MOI.supports(::NLToPolynomial, ::MOI.NLPBlock) = true
function MOI.set(model::NLToPolynomial{T}, ::MOI.NLPBlock, data::MOI.NLPBlockData) where {T}
    # FIXME if a non-NLP objective is set afterwards, it might overwrite.
    # but let's not complicate as it will be fixed by
    # https://github.com/jump-dev/MathOptInterface.jl/issues/846
    MOI.initialize(data.evaluator, [:ExprGraph])
    if data.has_objective
        d = Dict{MOI.VariableIndex, VarType}()
        poly = _to_polynomial!(d, T, MOI.objective_expr(data.evaluator))
        variable_map = collect(d)
        sort!(variable_map, by = x -> x[2])
        variables = [x[1] for x in variable_map]
        obj = FuncType{T}(poly, variables)
        MOI.set(model.model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    end
end

# The boilerplate part: passing everything to the inner `.model`

MOI.is_empty(model::NLToPolynomial) = MOI.is_empty(model.model)
MOI.empty!(model::NLToPolynomial) = MOI.empty!(model.model)

MOI.supports_incremental_interface(model::NLToPolynomial) = MOI.supports_incremental_interface(model.model)
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

function MOI.set(
    model::NLToPolynomial,
    attr::MOI.AbstractModelAttribute,
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
