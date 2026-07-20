"""
    abstract type AbstractPolynomialOptimizer{T} <: MOI.AbstractOptimizer end

Optimizer for polynomial optimization problems storing the problem in a
[`Model`](@ref).

Subtypes should be mutable structs with the fields
```julia
model::PolyJuMP.Model{T}
solve_time::Float64
```
and implement [`_invalidate!`](@ref) and [`_optimize!`](@ref).
"""
abstract type AbstractPolynomialOptimizer{T} <: MOI.AbstractOptimizer end

"""
    _invalidate!(model::AbstractPolynomialOptimizer)

Invalidate the result of previous calls to `MOI.optimize!`; called when the
polynomial optimization problem is modified.
"""
function _invalidate! end

"""
    _optimize!(model::AbstractPolynomialOptimizer)

Solve the polynomial optimization problem stored in `model.model`; called by
`MOI.optimize!` which records the elapsed time in `model.solve_time`.
"""
function _optimize! end

MOI.is_empty(model::AbstractPolynomialOptimizer) = MOI.is_empty(model.model)

function MOI.empty!(model::AbstractPolynomialOptimizer)
    MOI.empty!(model.model)
    _invalidate!(model)
    return
end

function MOI.supports(
    model::AbstractPolynomialOptimizer,
    attr::MOI.AbstractModelAttribute,
)
    return MOI.supports(model.model, attr)
end

function MOI.set(
    model::AbstractPolynomialOptimizer,
    attr::MOI.AbstractModelAttribute,
    value,
)
    MOI.set(model.model, attr, value)
    _invalidate!(model)
    return
end

function MOI.get(
    model::AbstractPolynomialOptimizer,
    attr::Union{
        MOI.AbstractModelAttribute,
        MOI.Bridges.ListOfNonstandardBridges,
    },
)
    return MOI.get(model.model, attr)
end

function MOI.is_valid(model::AbstractPolynomialOptimizer, i::MOI.Index)
    return MOI.is_valid(model.model, i)
end

function MOI.add_variable(model::AbstractPolynomialOptimizer)
    _invalidate!(model)
    return MOI.add_variable(model.model)
end

function MOI.supports_constraint(
    model::AbstractPolynomialOptimizer,
    ::Type{F},
    ::Type{S},
) where {F<:MOI.AbstractFunction,S<:MOI.AbstractSet}
    return MOI.supports_constraint(model.model, F, S)
end

function MOI.add_constraint(
    model::AbstractPolynomialOptimizer,
    func::MOI.AbstractFunction,
    set::MOI.AbstractSet,
)
    ci = MOI.add_constraint(model.model, func, set)
    _invalidate!(model)
    return ci
end

MOI.supports_incremental_interface(::AbstractPolynomialOptimizer) = true

function MOI.copy_to(dest::AbstractPolynomialOptimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

function MOI.optimize!(model::AbstractPolynomialOptimizer)
    model.solve_time = @elapsed _optimize!(model)
    return
end

function MOI.get(model::AbstractPolynomialOptimizer, ::MOI.SolveTimeSec)
    return model.solve_time
end

"""
    MultiplierMaxdegree()

A constraint attribute for the maximum degree of the multiplier of the
constraint in the certificate of nonnegativity of the Lagrangian used by an
[`AbstractRelaxationOptimizer`](@ref).
Increasing this degree gives a higher level of the hierarchy of relaxations,
hence a possibly tighter objective bound at the price of a larger relaxation.
By default, the degree `d - maxdegree(g)` is used for a constraint of
polynomial `g` where `d` is the smallest even number larger than the maximum
degree of the objective and constraint polynomials.
"""
struct MultiplierMaxdegree <: MOI.AbstractConstraintAttribute end

function MOI.Bridges.Constraint.invariant_under_function_conversion(
    ::MultiplierMaxdegree,
)
    return true
end

"""
    abstract type AbstractRelaxationOptimizer{T} <: AbstractPolynomialOptimizer{T} end

Optimizer computing a bound on the objective value of a polynomial
optimization problem
```
min  f(x)
s.t. g_i(x) ≥ 0
     h_j(x) = 0
```
by solving the relaxation
```
max  t
s.t. f - t - Σ_i σ_i * g_i - Σ_j μ_j * h_j ∈ C
     σ_i ∈ C
```
where `μ_j` are free polynomials and `C` is the cone of certified nonnegative
polynomials returned by [`nonnegativity_cone`](@ref)
(and conversely for a `max` problem). The bound is the objective value of this
relaxation and is returned as the `MOI.ObjectiveBound`; no primal solution is
computed so the `MOI.ResultCount` is zero. The degrees of the multipliers
`σ_i` and `μ_j` are given by the [`MultiplierMaxdegree`](@ref) constraint
attribute.

Subtypes should be mutable structs with the fields
```julia
model::PolyJuMP.Model{T}
multiplier_maxdegree::Dict{MOI.ConstraintIndex,Int}
solver::Any
relaxation::Union{Nothing,JuMP.GenericModel{T}}
solve_time::Float64
```
and implement [`nonnegativity_cone`](@ref).
"""
abstract type AbstractRelaxationOptimizer{T} <: AbstractPolynomialOptimizer{T} end

"""
    nonnegativity_cone(model::AbstractRelaxationOptimizer)

Return the set (e.g., `PolyJuMP.SAGE.Polynomials()` or
`SumOfSquares.SOSCone()`) in which a polynomial is constrained to belong as a
sufficient condition for its nonnegativity.
"""
function nonnegativity_cone end

function _invalidate!(model::AbstractRelaxationOptimizer)
    model.relaxation = nothing
    model.solve_time = NaN
    return
end

function MOI.empty!(model::AbstractRelaxationOptimizer)
    MOI.empty!(model.model)
    empty!(model.multiplier_maxdegree)
    _invalidate!(model)
    return
end

function MOI.supports(
    ::AbstractRelaxationOptimizer{T},
    ::MultiplierMaxdegree,
    ::Type{<:MOI.ConstraintIndex{<:ScalarPolynomialFunction{T}}},
) where {T}
    return true
end

function MOI.set(
    model::AbstractRelaxationOptimizer{T},
    ::MultiplierMaxdegree,
    ci::MOI.ConstraintIndex{<:ScalarPolynomialFunction{T}},
    degree::Integer,
) where {T}
    MOI.throw_if_not_valid(model, ci)
    model.multiplier_maxdegree[ci] = degree
    _invalidate!(model)
    return
end

function MOI.get(
    model::AbstractRelaxationOptimizer{T},
    ::MultiplierMaxdegree,
    ci::MOI.ConstraintIndex{<:ScalarPolynomialFunction{T}},
) where {T}
    return get(model.multiplier_maxdegree, ci, nothing)
end

_equalities(::SS.FullSpace) = []
_equalities(set::SS.AbstractAlgebraicSet) = SS.equalities(set)
_equalities(set::SS.BasicSemialgebraicSet) = _equalities(set.V)
_inequalities(::SS.AbstractAlgebraicSet) = []
_inequalities(set::SS.BasicSemialgebraicSet) = SS.inequalities(set)

function _multiplier(relaxation::JuMP.GenericModel, x, degree)
    basis = MB.SubBasis{MB.Monomial}(MP.monomials(x, 0:degree))
    poly = JuMP.@variable(relaxation, variable_type = Poly(basis))
    return MP.polynomial(poly)
end

function _optimize!(model::AbstractRelaxationOptimizer{T}) where {T}
    pop = model.model
    x = MP.variables(pop)
    if pop.objective_sense == MOI.FEASIBILITY_SENSE ||
       isnothing(pop.objective_function)
        f = zero(PolyType{T})
    else
        f = pop.objective_function
    end
    eqs = _equalities(pop.set)
    ineqs = _inequalities(pop.set)
    maxdeg = max(
        MP.maxdegree(f),
        maximum(MP.maxdegree, eqs; init = 0),
        maximum(MP.maxdegree, ineqs; init = 0),
    )
    maxdeg += isodd(maxdeg)
    eq_degrees = Dict{Int,Int}()
    ineq_degrees = Dict{Int,Int}()
    for (ci, d) in model.multiplier_maxdegree
        if ci isa MOI.ConstraintIndex{<:Any,<:MOI.EqualTo}
            eq_degrees[ci.value] = d
        else
            ineq_degrees[ci.value] = d
        end
    end
    relaxation = JuMP.GenericModel{T}(model.solver)
    cone = nonnegativity_cone(model)
    t = JuMP.@variable(relaxation, base_name = "t")
    lagrangian = f - t
    if pop.objective_sense == MOI.MAX_SENSE
        lagrangian = -lagrangian
    end
    for (i, h) in enumerate(eqs)
        d = get(eq_degrees, i, maxdeg - MP.maxdegree(h))
        lagrangian -= _multiplier(relaxation, x, d) * h
    end
    for (i, g) in enumerate(ineqs)
        d = get(ineq_degrees, i, maxdeg - MP.maxdegree(g))
        σ = _multiplier(relaxation, x, d)
        JuMP.@constraint(relaxation, σ in cone)
        lagrangian -= σ * g
    end
    JuMP.@constraint(relaxation, lagrangian in cone)
    sense = pop.objective_sense == MOI.MAX_SENSE ? MOI.MIN_SENSE : MOI.MAX_SENSE
    JuMP.set_objective(relaxation, sense, t)
    JuMP.optimize!(relaxation)
    model.relaxation = relaxation
    return
end

function MOI.get(model::AbstractRelaxationOptimizer, ::MOI.TerminationStatus)
    if isnothing(model.relaxation)
        return MOI.OPTIMIZE_NOT_CALLED
    end
    return JuMP.termination_status(model.relaxation)
end

function MOI.get(model::AbstractRelaxationOptimizer, ::MOI.RawStatusString)
    if isnothing(model.relaxation)
        return "`optimize!` has not yet been called"
    end
    return JuMP.raw_status(model.relaxation)
end

function MOI.get(model::AbstractRelaxationOptimizer, ::MOI.ObjectiveBound)
    return JuMP.objective_value(model.relaxation)
end

MOI.get(::AbstractRelaxationOptimizer, ::MOI.ResultCount) = 0

MOI.get(::AbstractRelaxationOptimizer, ::MOI.PrimalStatus) = MOI.NO_SOLUTION

MOI.get(::AbstractRelaxationOptimizer, ::MOI.DualStatus) = MOI.NO_SOLUTION
