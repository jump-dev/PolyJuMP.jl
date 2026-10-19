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

# The matrix inequalities stored in the `Model` are only supported by the
# optimizers implementing `matrix_nonnegativity_cone`.
function MOI.supports_constraint(
    model::AbstractPolynomialOptimizer,
    ::Type{F},
    ::Type{S},
) where {
    F<:MOI.AbstractVectorFunction,
    S<:MOI.AbstractSymmetricMatrixSetTriangle,
}
    return !isnothing(matrix_nonnegativity_cone(model)) &&
           MOI.supports_constraint(model.model, F, S)
end

# The bridging costs are forwarded to `model.model` with the other model
# attributes but they should be consistent with the constraints supported by
# `model`.
function MOI.get(
    model::AbstractPolynomialOptimizer,
    ::MOI.ConstraintBridgingCost{F,S},
) where {
    F<:MOI.AbstractVectorFunction,
    S<:MOI.AbstractSymmetricMatrixSetTriangle,
}
    return MOI.supports_constraint(model, F, S) ? 0.0 : Inf
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
     G_k(x) ⪰ 0
```
by solving the relaxation
```
max  t
s.t. f - t - Σ_i σ_i * g_i - Σ_j μ_j * h_j - Σ_k ⟨Σ_k, G_k⟩ ∈ C
     σ_i ∈ C
     Σ_k ∈ M
```
where `μ_j` are free polynomials, `C` is the cone of certified nonnegative
polynomials returned by [`nonnegativity_cone`](@ref) and `M` is the cone of
certified positive semidefinite polynomial matrices returned by
[`matrix_nonnegativity_cone`](@ref)
(and conversely for a `max` problem). The bound is the objective value of this
relaxation and is returned as the `MOI.ObjectiveBound`. The degrees of the
multipliers `σ_i`, `μ_j` and of the entries of `Σ_k` are given by the
[`MultiplierMaxdegree`](@ref) constraint attribute.
The polynomial matrix inequalities `G_k(x) ⪰ 0` are only supported if
[`matrix_nonnegativity_cone`](@ref) is implemented, see [HL06] for the
Sum-of-Squares case. Candidate primal solutions may in addition be recovered
from the solution of the relaxation by implementing
[`recover_solutions`](@ref); the `MOI.ResultCount` is the number of
recovered candidates.

Subtypes should be mutable structs with the fields
```julia
model::PolyJuMP.Model{T}
multiplier_maxdegree::Dict{MOI.ConstraintIndex,Int}
solver::Any
relaxation::Union{Nothing,JuMP.GenericModel{T}}
solutions::Vector{PolyJuMP.Solution{T}}
feasibility_tolerance::T
solve_time::Float64
```
and implement [`nonnegativity_cone`](@ref) and [`recover_solutions`](@ref).

[HL06] Henrion, Didier and Lasserre, Jean-Bernard.
*Convergent relaxations of polynomial matrix inequalities and static output
feedback.*
IEEE Transactions on Automatic Control 51.2 (2006): 192-202.
"""
abstract type AbstractRelaxationOptimizer{T} <: AbstractPolynomialOptimizer{T} end

"""
    nonnegativity_cone(model::AbstractRelaxationOptimizer)

Return the set (e.g., `PolyJuMP.SAGE.Polynomials()` or
`SumOfSquares.SOSCone()`) in which a polynomial is constrained to belong as a
sufficient condition for its nonnegativity.
"""
function nonnegativity_cone end

"""
    matrix_nonnegativity_cone(model::AbstractPolynomialOptimizer)

Return the set (e.g., `SumOfSquares.SOSMatrixCone()`) in which a symmetric
matrix of polynomials is constrained to belong as a sufficient condition for
its positive semidefiniteness, or `nothing` (the default) if polynomial matrix
inequalities are not supported by `model`.
"""
matrix_nonnegativity_cone(::AbstractPolynomialOptimizer) = nothing

function _invalidate!(model::AbstractRelaxationOptimizer)
    model.relaxation = nothing
    empty!(model.solutions)
    model.solve_time = NaN
    return
end

"""
    recover_solutions(
        model::AbstractRelaxationOptimizer{T},
        relaxation::JuMP.GenericModel{T},
        cref::JuMP.ConstraintRef,
        lagrangian,
    ) where {T}

Return a vector of candidate [`Solution`](@ref)s recovered from the solution
of `relaxation`, where `cref` is the constraint of the `lagrangian` polynomial
in the cone [`nonnegativity_cone`](@ref). Return an empty vector if no
solution is recovered; `PolyJuMP.SAGE.Optimizer` implements the recovery from
the dual of `cref` of [MCW21, Section 4.2].

[MCW21] Murray, Riley, Venkat Chandrasekaran, and Adam Wierman.
"Signomials and polynomial optimization via relative entropy and partial dualization."
Mathematical Programming Computation 13 (2021): 257-295.
https://arxiv.org/pdf/1907.00814.pdf
"""
function recover_solutions end

function MOI.empty!(model::AbstractRelaxationOptimizer)
    MOI.empty!(model.model)
    empty!(model.multiplier_maxdegree)
    _invalidate!(model)
    return
end

const _MatrixIndex = MOI.ConstraintIndex{
    <:MOI.AbstractVectorFunction,
    <:MOI.AbstractSymmetricMatrixSetTriangle,
}

const _MultiplierIndex{T} =
    Union{MOI.ConstraintIndex{<:ScalarPolynomialFunction{T}},_MatrixIndex}

function MOI.supports(
    ::AbstractRelaxationOptimizer{T},
    ::MultiplierMaxdegree,
    ::Type{<:_MultiplierIndex{T}},
) where {T}
    return true
end

function MOI.set(
    model::AbstractRelaxationOptimizer{T},
    ::MultiplierMaxdegree,
    ci::_MultiplierIndex{T},
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
    ci::_MultiplierIndex{T},
) where {T}
    return get(model.multiplier_maxdegree, ci, nothing)
end

# The matrix inequalities in a square set are bridged into a triangle set so
# their multiplier is the multiplier of the triangle constraint.
function MOI.supports(
    model::MOI.ModelLike,
    attr::MultiplierMaxdegree,
    ::Type{<:MOI.Bridges.Constraint.SquareBridge{T,F,G,TT}},
) where {T,F,G,TT}
    return MOI.supports(model, attr, MOI.ConstraintIndex{F,TT})
end

function MOI.set(
    model::MOI.ModelLike,
    attr::MultiplierMaxdegree,
    bridge::MOI.Bridges.Constraint.SquareBridge,
    degree,
)
    MOI.set(model, attr, bridge.triangle, degree)
    return
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MultiplierMaxdegree,
    bridge::MOI.Bridges.Constraint.SquareBridge,
)
    return MOI.get(model, attr, bridge.triangle)
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

_maxdegree(G::Matrix) = maximum(MP.maxdegree, G)

function _matrix_multiplier(relaxation::JuMP.GenericModel, x, n, degree)
    Σ = Matrix{MP.polynomial_type(VarType, JuMP.AffExpr)}(undef, n, n)
    for j in 1:n, i in 1:j
        Σ[i, j] = Σ[j, i] = _multiplier(relaxation, x, degree)
    end
    return Σ
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
    mats = pop.matrix_inequalities
    maxdeg = max(
        MP.maxdegree(f),
        maximum(MP.maxdegree, eqs; init = 0),
        maximum(MP.maxdegree, ineqs; init = 0),
        maximum(_maxdegree, mats; init = 0),
    )
    maxdeg += isodd(maxdeg)
    eq_degrees = Dict{Int,Int}()
    ineq_degrees = Dict{Int,Int}()
    mat_degrees = Dict{Int,Int}()
    for (ci, d) in model.multiplier_maxdegree
        if ci isa _MatrixIndex
            mat_degrees[ci.value] = d
        elseif ci isa MOI.ConstraintIndex{<:Any,<:MOI.EqualTo}
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
    for (i, G) in enumerate(mats)
        d = get(mat_degrees, i, maxdeg - _maxdegree(G))
        Σ = _matrix_multiplier(relaxation, x, size(G, 1), d)
        JuMP.@constraint(relaxation, Σ in matrix_nonnegativity_cone(model))
        lagrangian -= sum(G .* Σ)
    end
    cref = JuMP.@constraint(relaxation, lagrangian in cone)
    sense = pop.objective_sense == MOI.MAX_SENSE ? MOI.MIN_SENSE : MOI.MAX_SENSE
    JuMP.set_objective(relaxation, sense, t)
    JuMP.optimize!(relaxation)
    model.relaxation = relaxation
    model.solutions = recover_solutions(model, relaxation, cref, lagrangian)
    postprocess!(model.solutions, pop, nothing)
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

function MOI.get(model::AbstractRelaxationOptimizer, ::MOI.ResultCount)
    return length(model.solutions)
end

function MOI.get(model::AbstractRelaxationOptimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return model.solutions[attr.result_index].objective_value
end

function MOI.get(
    model::AbstractRelaxationOptimizer,
    attr::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.throw_if_not_valid(model, vi)
    MOI.check_result_index_bounds(model, attr)
    return model.solutions[attr.result_index].values[vi.value]
end

function MOI.get(model::AbstractRelaxationOptimizer, attr::MOI.PrimalStatus)
    if attr.result_index in 1:length(model.solutions)
        return model.solutions[attr.result_index].status
    end
    return MOI.NO_SOLUTION
end

MOI.get(::AbstractRelaxationOptimizer, ::MOI.DualStatus) = MOI.NO_SOLUTION
