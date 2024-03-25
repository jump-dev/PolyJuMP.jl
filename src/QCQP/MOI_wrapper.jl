import MathOptInterface as MOI

const _Model{F,S} =
    MOI.Utilities.UniversalFallback{MOI.Utilities.VectorOfConstraints{F,S}}

mutable struct Optimizer{T,O<:MOI.ModelLike} <: MOI.AbstractOptimizer
    model::O
    objective::Union{Nothing,PolyJuMP.ScalarPolynomialFunction{T}}
    constraints::DataStructures.OrderedDict{Type,Tuple{Type,_Model}}
    index_map::Dict{MOI.ConstraintIndex,MOI.ConstraintIndex}
end

function Optimizer{T}(model::MOI.ModelLike) where {T}
    return Optimizer{T,typeof(model)}(
        model,
        nothing,
        DataStructures.OrderedDict{Type,MOI.Utilities.VectorOfConstraints}(),
        Dict{MOI.ConstraintIndex,MOI.ConstraintIndex}(),
    )
end

Optimizer(model::MOI.ModelLike) = Optimizer{Float64}(model)

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
function MOI.empty!(model::Optimizer)
    MOI.empty!(model.model)
    model.objective = nothing
    empty!(model.constraints)
    empty!(model.index_map)
    return
end

MOI.is_valid(model::Optimizer, i::MOI.Index) = MOI.is_valid(model.model, i)
function MOI.is_valid(
    model::Optimizer{T},
    ::MOI.ConstraintIndex{PolyJuMP.ScalarPolynomialFunction{T},S},
) where {T,S<:MOI.AbstractScalarSet}
    return haskey(model.constraints, S) &&
           MOI.is_valid(model.constraints[S][2], ci)
end

function MOI.supports(
    model::Optimizer,
    attr::MOI.AbstractConstraintAttribute,
    C::Type{<:MOI.ConstraintIndex},
)
    return MOI.supports(model.model, attr, C)
end

function MOI.set(
    model::Optimizer,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
    value,
)
    return MOI.set(model.model, attr, ci, value)
end

function MOI.set(
    model::Optimizer{T},
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex{<:PolyJuMP.ScalarPolynomialFunction{T},S},
    value,
) where {T,S<:MOI.AbstractScalarSet}
    return MOI.get(model.constraints[S][2], attr, model.index_map[ci], value)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
)
    return MOI.get(model.model, attr, ci)
end

function MOI.get(
    model::Optimizer{T},
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex{<:PolyJuMP.ScalarPolynomialFunction{T},S},
) where {T,S<:MOI.AbstractScalarSet}
    if MOI.is_set_by_optimize(attr)
        return MOI.get(model.model, attr, model.index_map[ci])
    else
        return MOI.get(model.constraints[S][2], attr, ci)
    end
end

MOI.add_variable(model::Optimizer) = MOI.add_variable(model.model)

function MOI.supports(
    model::Optimizer,
    attr::MOI.AbstractVariableAttribute,
    ::Type{MOI.VariableIndex},
)
    return MOI.supports(model.model, attr, MOI.VariableIndex)
end

function MOI.set(
    model::Optimizer,
    attr::MOI.AbstractVariableAttribute,
    vi::MOI.VariableIndex,
    value,
)
    return MOI.set(model.model, attr, vi, value)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.AbstractVariableAttribute,
    vi::MOI.VariableIndex,
)
    return MOI.get(model.model, attr, vi)
end

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
    return MOI.set(model.model, attr, value)
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

function MOI.supports_constraint(
    model::Optimizer{T},
    ::Type{<:PolyJuMP.ScalarPolynomialFunction{T}},
    ::Type{S},
) where {T,S<:MOI.AbstractScalarSet}
    return MOI.supports_constraint(
        model.model,
        MOI.ScalarQuadraticFunction{T},
        S,
    )
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.AbstractFunction,
    set::MOI.AbstractSet,
)
    return MOI.add_constraint(model.model, func, set)
end

function MOI.add_constraint(
    model::Optimizer{T},
    func::PolyJuMP.ScalarPolynomialFunction{T,P},
    set::MOI.AbstractScalarSet,
) where {T,P}
    F = typeof(func)
    S = typeof(set)
    if !haskey(model.constraints, S)
        con = MOI.Utilities.UniversalFallback(
            MOI.Utilities.VectorOfConstraints{F,S}(),
        )
        model.constraints[S] = (P, con)
    end
    return MOI.add_constraint(model.constraints[S][2], func, set)
end

function MOI.get(
    model::Optimizer{T},
    attr::Union{MOI.ConstraintFunction,MOI.ConstraintSet},
    ci::MOI.ConstraintIndex{<:PolyJuMP.ScalarPolynomialFunction{T},S},
) where {T,S<:MOI.AbstractScalarSet}
    return MOI.get(model.constraints[S][2], attr, ci)
end

function MOI.get(
    model::Optimizer{T},
    attr::MOI.ListOfConstraintIndices{<:PolyJuMP.ScalarPolynomialFunction{T},S},
) where {T,S<:MOI.AbstractScalarSet}
    return MOI.get(model.constraints[S][2], attr)
end

function MOI.supports_incremental_interface(model::Optimizer)
    return MOI.supports_incremental_interface(model.model)
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

MOI.optimize!(model::Optimizer) = MOI.optimize!(model.model)

function _quad_convert(p::MP.AbstractPolynomialLike{T}, index, div) where {T}
    q = zero(MOI.ScalarQuadraticFunction{T})
    for t in MP.terms(p)
        α = MP.coefficient(t)
        mono = MP.monomial(t)
        if MP.degree(mono) == 0
            MA.operate!(+, q, α)
        else
            if haskey(index, mono)
                MA.operate!(MA.add_mul, q, α, index[mono])
            else
                x = div[mono]
                y = MP.div_multiple(mono, x)
                MA.operate!(MA.add_mul, q, α, index[x], index[y])
            end
        end
    end
    return q
end

function _add_monomials!(p::PolyJuMP.ScalarPolynomialFunction, monos1)
    monos2 = MP.monomials(p.polynomial)
    if isnothing(monos1)
        return monos2
    else
        return MP.merge_monomial_vectors([monos1, monos2])
    end
end

function _subs!(
    p::PolyJuMP.ScalarPolynomialFunction{T,P},
    ::Nothing,
) where {T,P}
    return p,
    Dict{MOI.VariableIndex,MP.variable_union_type(P)}(
        vi => var for (vi, var) in zip(p.variables, MP.variables(p.polynomial))
    )
end

"""
    _subs_ensure_moi_order(p::PolyJuMP.ScalarPolynomialFunction, old, new)

Substitutes old `MP.variables(p.polynomial)` with new vars, while re-sorting the
MOI `p.variables` to get them in the correct order after substitution.
"""
function _subs_ensure_moi_order(p::PolyJuMP.ScalarPolynomialFunction, old, new)
    if isempty(old)
        return p
    end
    poly = MP.subs(p.polynomial, old => new)
    all_new_vars = MP.variables(poly)
    to_old_map = Dict(zip(new, old))
    to_moi_map = Dict(zip(MP.variables(p.polynomial), p.variables))
    moi_vars = [to_moi_map[get(to_old_map, v, v)] for v in all_new_vars]
    return PolyJuMP.ScalarPolynomialFunction(poly, moi_vars)
end

function _subs!(
    p::PolyJuMP.ScalarPolynomialFunction,
    index_to_var::Dict{K,V},
) where {K,V}
    old_var = V[]
    new_var = V[]
    for (vi, var) in zip(p.variables, MP.variables(p.polynomial))
        if haskey(index_to_var, vi)
            if var != index_to_var[vi]
                push!(old_var, var)
                push!(new_var, index_to_var[vi])
            end
        else
            index_to_var[vi] = var
        end
    end
    p = _subs_ensure_moi_order(p, old_var, new_var)
    return p, index_to_var
end

function _add_variables!(
    p::PolyJuMP.ScalarPolynomialFunction{T,P},
    d,
) where {T,P}
    if isnothing(d)
        d = Dict{MP.monomial_type(P),MOI.VariableIndex}()
    else
        M = promote_type(keytype(d), MP.monomial_type(P))
        if keytype(d) !== M
            d = convert(Dict{M,MOI.VariableIndex}, d)
        end
    end
    for (v, vi) in zip(MP.variables(p.polynomial), p.variables)
        d[v] = vi
    end
    return d
end

import IntervalArithmetic

function monomial_variable_index(
    model::Optimizer{T},
    d::Dict,
    div,
    mono::MP.AbstractMonomialLike,
) where {T}
    if !haskey(d, mono)
        # If we don't have a variable for `mono` yet,
        # we create one now by equal to `x * y`.
        mono_bounds = IntervalArithmetic.interval(one(T))
        mono_start = one(T)
        for var in MP.variables(mono)
            deg = MP.degree(mono, var)
            if deg == 0
                continue
            end
            vi = monomial_variable_index(model, d, div, var)
            lb_var, ub_var = MOI.Utilities.get_bounds(model, T, vi)
            F = float(T)
            var_bounds = IntervalArithmetic.interval(
                lb_var == typemin(lb_var) ? typemin(F) : float(lb_var),
                ub_var == typemax(ub_var) ? typemax(F) : float(ub_var),
            )
            mono_bounds *= var_bounds^deg
            attr = MOI.VariablePrimalStart()
            if !isnothing(mono_start) &&
               MOI.supports(model, attr, MOI.VariableIndex)
                start = MOI.get(model, attr, vi)
                if isnothing(start)
                    mono_start = nothing
                else
                    mono_start *= start^deg
                end
            end
        end
        x = div[mono]
        vx = monomial_variable_index(model, d, div, x)
        y = MP.div_multiple(mono, x)
        vy = monomial_variable_index(model, d, div, y)
        lb = IntervalArithmetic.inf(mono_bounds)
        ub = IntervalArithmetic.sup(mono_bounds)
        if isfinite(lb)
            if isfinite(ub)
                if lb == ub
                    set = MOI.EqualTo(T(lb))
                else
                    set = MOI.Interval(T(lb), T(ub))
                end
            else
                set = MOI.GreaterThan(T(lb))
            end
        else
            if isfinite(ub)
                set = MOI.LessThan(T(ub))
            else
                set = nothing
            end
        end
        if isnothing(set)
            d[mono] = MOI.add_variable(model.model)
        else
            d[mono], _ = MOI.add_constrained_variable(model.model, set)
        end
        if !isnothing(mono_start) &&
           MOI.supports(model, MOI.VariablePrimalStart(), MOI.VariableIndex)
            MOI.set(model.model, MOI.VariablePrimalStart(), d[mono], mono_start)
        end
        MOI.Utilities.normalize_and_add_constraint(
            model,
            MA.@rewrite(one(T) * d[mono] - one(T) * vx * vy),
            MOI.EqualTo(zero(T));
            allow_modify_function = true,
        )
    end
    return d[mono]
end

function _add_constraints(
    dest,
    src,
    index_map,
    cis_src::Vector{MOI.ConstraintIndex{F,S}},
    index_to_var,
    d,
    div,
) where {F,S}
    for ci in cis_src
        func = MOI.get(src, MOI.ConstraintFunction(), ci)
        set = MOI.get(src, MOI.ConstraintSet(), ci)
        func, index_to_var = _subs!(func, index_to_var)
        quad = _quad_convert(func.polynomial, d, div)
        dest_ci = MOI.Utilities.normalize_and_add_constraint(dest, quad, set)
        index_map[ci] = dest_ci
    end
    # `Utilities.pass_attributes` needs `index_map` to be an `IndexMap` :(
    #MOI.Utilities.pass_attributes(dest, src, index_map, cis_src)
    # `ListOfConstraintAttributesSet` not defined for `VectorOfConstraints`
    #    for attr in MOI.get(src, MOI.ListOfConstraintAttributesSet{F,S}())
    #        if !MOI.supports(dest, attr)
    #            if attr == MOI.Name()
    #                continue  # Skipping names is okay.
    #            end
    #        end
    #        for ci in cis_src
    #            value = MOI.get(src, attr, ci)
    #            if value !== nothing
    #                MOI.set(
    #                    dest,
    #                    attr,
    #                    index_map[ci],
    #                    MOI.Utilities.map_indices(index_map, attr, value),
    #                )
    #            end
    #        end
    #    end
    return
end

function MOI.Utilities.final_touch(model::Optimizer{T}, _) where {T}
    index_to_var = nothing
    vars = nothing
    monos = nothing
    if !isnothing(model.objective)
        func, index_to_var = _subs!(model.objective, index_to_var)
        vars = _add_variables!(func, vars)
        monos = _add_monomials!(func, monos)
    end
    for (S, constraints) in model.constraints
        F = PolyJuMP.ScalarPolynomialFunction{T,constraints[1]}
        for ci in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
            func = MOI.get(model, MOI.ConstraintFunction(), ci)
            func, index_to_var = _subs!(func, index_to_var)
            vars = _add_variables!(func, vars)
            monos = _add_monomials!(func, monos)
        end
    end
    if !isnothing(monos)
        div = decompose(monos)
        for mono in sort(collect(keys(div)))
            if haskey(vars, mono)
                continue
            end
            a = div[mono]
            monomial_variable_index(model, vars, div, a)
            b = MP.div_multiple(mono, a)
            monomial_variable_index(model, vars, div, b)
        end
        if !isnothing(model.objective)
            func, index_to_var = _subs!(model.objective, index_to_var)
            obj = _quad_convert(func.polynomial, vars, div)
            MOI.set(model.model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
        end
        for S in keys(model.constraints)
            F = PolyJuMP.ScalarPolynomialFunction{T,model.constraints[S][1]}
            cis = MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
            src = model.constraints[S][2]
            _add_constraints(
                model.model,
                src,
                model.index_map,
                cis,
                index_to_var,
                vars,
                div,
            )
        end
    end
    return
end

function MOI.get(model::Optimizer, attr::MOI.SolverName)
    name = MOI.get(model.model, attr)
    return "PolyJuMP.QCQP with $name"
end
