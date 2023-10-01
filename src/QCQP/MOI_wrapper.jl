import MathOptInterface as MOI

mutable struct Optimizer{T,O<:MOI.ModelLike} <: MOI.AbstractOptimizer
    model::O
    objective::Union{Nothing,PolyJuMP.ScalarPolynomialFunction{T}}
    constraints::DataStructures.OrderedDict{Type,Tuple{Type,MOI.Utilities.VectorOfConstraints}}
end

function Optimizer{T}(model::MOI.ModelLike) where {T}
    return Optimizer{T,typeof(model)}(
        model,
        nothing,
        DataStructures.OrderedDict{Type,MOI.Utilities.VectorOfConstraints}(),
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
    return
end

MOI.is_valid(model::Optimizer, i::MOI.Index) = MOI.is_valid(model.model, i)
function MOI.is_valid(
    model::Optimizer{T},
    ::MOI.ConstraintIndex{
        PolyJuMP.ScalarPolynomialFunction{T},
        S
    },
) where {T,S<:MOI.AbstractScalarSet}
    return haskey(model.constraints, S) &&
        MOI.is_valid(model.constraints[S][2], ci)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
)
    return MOI.get(model.model, attr, ci)
end

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
        con = MOI.Utilities.VectorOfConstraints{F,S}()
        model.constraints[S] = (P,con)
    end
    return MOI.add_constraint(model.constraints[S][2], func, set)
end

function MOI.get(
    model::Optimizer{T},
    attr::Union{
        MOI.ConstraintFunction,
        MOI.ConstraintSet,
    },
    ci::MOI.ConstraintIndex{
        <:PolyJuMP.ScalarPolynomialFunction{T},
        S,
    },
) where {T,S}
    return MOI.get(model.constraints[S][2], attr, ci)
end

function MOI.get(
    model::Optimizer{T},
    attr::MOI.ListOfConstraintIndices{
        <:PolyJuMP.ScalarPolynomialFunction{T},
        S,
    }
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

function _subs!(p::PolyJuMP.ScalarPolynomialFunction{T,P}, ::Nothing) where {T,P}
    return p, Dict{MOI.VariableIndex,MP.variable_union_type(P)}(
        vi => var for (vi, var) in zip(p.variables, MP.variables(p.polynomial))
    )
end

function _subs!(p::PolyJuMP.ScalarPolynomialFunction, index_to_var::Dict{K,V}) where {K,V}
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
    if !isempty(old_var)
        poly = MP.subs(p.polynomial, old_var => new_var)
        p = PolyJuMP.ScalarPolynomialFunction(poly, p.variables)
    end
    return p, index_to_var
end

function _add_variables!(p::PolyJuMP.ScalarPolynomialFunction{T,P}, d) where {T,P}
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

function monomial_variable_index(model::Optimizer{T}, d::Dict, div, mono::MP.AbstractMonomialLike) where {T}
    if !haskey(d, mono)
        x = div[mono]
        vx = monomial_variable_index(model, d, div, x)
        y = MP.div_multiple(mono, x)
        vy = monomial_variable_index(model, d, div, y)
        lx, ux = MOI.Utilities.get_bounds(model, T, vx)
        ly, uy = MOI.Utilities.get_bounds(model, T, vy)
        bounds = (lx * ly, lx * uy, ux * ly, ux * uy)
        l = min(bounds...)
        if vx == vy
            l = max(l, zero(T))
        end
        u = max(bounds...)
        d[mono], _ = MOI.add_constrained_variable(model.model, MOI.Interval(l, u))
        MOI.add_constraint(model,
            MA.@rewrite(one(T) * d[mono] - one(T) * vx * vy),
            MOI.EqualTo(zero(T)),
        )
    end
    return d[mono]
end

function _add_constraints(model::Optimizer, cis, index_to_var, d, div)
    for ci in cis
        func = MOI.get(model, MOI.ConstraintFunction(), ci)
        set = MOI.get(model, MOI.ConstraintSet(), ci)
        func, index_to_var = _subs!(func, index_to_var)
        quad = _quad_convert(func.polynomial, d, div)
        MOI.add_constraint(model, quad, set)
    end
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
    if !isempty(model.constraints)
        for S in keys(model.constraints)
            for ci in MOI.get(model, MOI.ListOfConstraintIndices{
                PolyJuMP.ScalarPolynomialFunction{T,model.constraints[S][1]},
                S
            }())
                func = MOI.get(model, MOI.ConstraintFunction(), ci)
                func, index_to_var = _subs!(func, index_to_var)
                vars = _add_variables!(func, vars)
                monos = _add_monomials!(func, monos)
            end
        end
    end
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
        _add_constraints(model, cis, index_to_var, vars, div)
    end
    return
end

function MOI.get(model::Optimizer, attr::MOI.SolverName)
    name = MOI.get(model.model, attr)
    return "PolyJuMP.QCQP with $name"
end