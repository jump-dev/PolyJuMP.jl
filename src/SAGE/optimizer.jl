"""
    Optimizer{T}(solver; feasibility_tolerance = sqrt(Base.rtoldefault(T)))

Optimizer computing a bound on the objective value of a polynomial
optimization problem using its SAGE relaxation, see
[`PolyJuMP.AbstractRelaxationOptimizer`](@ref) with the cone
[`Polynomials`](@ref) as certificate of nonnegativity [CP16; MCW21].
The relaxation is solved with `solver` and the bound can be queried with
`MOI.ObjectiveBound()`. Candidate solutions are in addition recovered from
the dual of the relaxation following [MCW21, Section 4.2]; the candidates
are classified as feasible up to `feasibility_tolerance` (see
`MOI.PrimalStatus`) and sorted by objective value.

[CP16] Chandrasekaran, Venkat, and Parikshit Shah.
"Relative entropy relaxations for signomial optimization."
SIAM Journal on Optimization 26.2 (2016): 1147-1173.
[MCW21] Murray, Riley, Venkat Chandrasekaran, and Adam Wierman.
"Signomials and polynomial optimization via relative entropy and partial dualization."
Mathematical Programming Computation 13 (2021): 257-295.
https://arxiv.org/pdf/1907.00814.pdf
"""
mutable struct Optimizer{T} <: PolyJuMP.AbstractRelaxationOptimizer{T}
    model::PolyJuMP.Model{T}
    multiplier_maxdegree::Dict{MOI.ConstraintIndex,Int}
    solver::Any
    relaxation::Union{Nothing,JuMP.GenericModel{T}}
    solutions::Vector{PolyJuMP.Solution{T}}
    feasibility_tolerance::T
    solve_time::Float64
end

function Optimizer{T}(
    solver;
    feasibility_tolerance = sqrt(Base.rtoldefault(T)),
) where {T}
    return Optimizer{T}(
        PolyJuMP.Model{T}(),
        Dict{MOI.ConstraintIndex,Int}(),
        solver,
        nothing,
        PolyJuMP.Solution{T}[],
        feasibility_tolerance,
        NaN,
    )
end

Optimizer(solver; kws...) = Optimizer{Float64}(solver; kws...)

MOI.get(::Optimizer, ::MOI.SolverName) = "PolyJuMP.SAGE"

PolyJuMP.nonnegativity_cone(::Optimizer) = Polynomials()

# Solve `A * z = b` over GF(2). Return `nothing` if the system is infeasible
# and otherwise a particular solution and a basis of the nullspace of `A`.
# This is used for the sign recovery of [MCW21, Section 4.2], see also
# `sageopt.relaxations.symbolic_correspondences.mod2linsolve` in `sageopt`:
# https://github.com/rileyjmurray/sageopt
function _mod2_solve(A::Matrix{Bool}, b::Vector{Bool})
    A = copy(A)
    b = copy(b)
    m, n = size(A)
    pivot_cols = Int[]
    for col in 1:n
        r = length(pivot_cols) + 1
        p = findfirst(i -> A[i, col], r:m)
        if isnothing(p)
            continue
        end
        p += r - 1
        A[[r, p], :] = A[[p, r], :]
        b[r], b[p] = b[p], b[r]
        for i in 1:m
            if i != r && A[i, col]
                A[i, :] .⊻= A[r, :]
                b[i] ⊻= b[r]
            end
        end
        push!(pivot_cols, col)
    end
    if any(i -> b[i], (length(pivot_cols)+1):m)
        return nothing
    end
    z = fill(false, n)
    for (i, col) in enumerate(pivot_cols)
        z[col] = b[i]
    end
    nullspace = map(setdiff(1:n, pivot_cols)) do col
        w = fill(false, n)
        w[col] = true
        for (i, pivot) in enumerate(pivot_cols)
            w[pivot] = A[i, col]
        end
        return w
    end
    return z, nullspace
end

"""
    PolyJuMP.recover_solutions(model::Optimizer, relaxation, cref, lagrangian)

Recover candidate solutions from the dual `v` of the SAGE constraint `cref`
of the Lagrangian, following [MCW21, Section 4.2] whose reference
implementation is `poly_solrec` in `sageopt`:
https://github.com/rileyjmurray/sageopt/blob/master/sageopt/relaxations/poly_solution_recovery.py

After normalizing `v` by its entry for the constant monomial, `v` is a
pseudo-moment vector: at the optimum of a tight relaxation, `v[i]` is the
monomial `monos[i]` evaluated at an optimal solution (or a convex combination
of the moments of several optimal solutions). The magnitude of the variables
is recovered by solving the least squares problem `α * y ≈ log.(abs.(v))` on
the exponents `α` of the monomials with nonnegligible moment; variables
appearing in no such monomial have all their moments negligible, hence
magnitude zero. The signs are recovered from the linear system
`α * z ≡ [v[i] < 0] (mod 2)` over GF(2); when the moments leave signs
undetermined (e.g., for sign-symmetric problems), one candidate per element
of the affine solution set is returned.
"""
function PolyJuMP.recover_solutions(
    model::Optimizer{T},
    relaxation::JuMP.GenericModel{T},
    cref::JuMP.ConstraintRef,
    lagrangian,
) where {T}
    solutions = PolyJuMP.Solution{T}[]
    if JuMP.termination_status(relaxation) != MOI.OPTIMAL
        return solutions
    end
    v = MOI.get(
        JuMP.backend(relaxation),
        MOI.ConstraintDual(),
        JuMP.index(cref),
    )
    monos = MP.monomials(lagrangian)
    # The constant monomial is always present since `lagrangian` contains `t`
    i0 = findfirst(iszero ∘ MP.degree, monos)
    v /= v[i0]
    ztol = sqrt(Base.rtoldefault(T))
    vars = MP.variables(lagrangian)
    rows = [i for i in eachindex(v) if i != i0 && abs(v[i]) > ztol]
    cols = filter(eachindex(vars)) do j
        return any(i -> !iszero(MP.degree(monos[i], vars[j])), rows)
    end
    F = float(T)
    A = F[MP.degree(monos[i], vars[j]) for i in rows, j in cols]
    b = F[log(abs(v[i])) for i in rows]
    y = A \ b
    mags = zeros(T, length(vars))
    for (k, j) in enumerate(cols)
        mags[j] = exp(y[k])
    end
    As = Bool[
        isodd(MP.degree(monos[i], vars[j])) for i in rows, j in eachindex(vars)
    ]
    bs = Bool[v[i] < 0 for i in rows]
    signs = _mod2_solve(As, bs)
    if isnothing(signs)
        # No consistent signs; only attempt nonnegative variable values,
        # similar to the `heuristic_signs` option of `poly_solrec` in `sageopt`
        patterns = [fill(false, length(vars))]
    else
        z, nullspace = signs
        # Flipping the sign of a variable of magnitude zero gives the same
        # solution so we do not enumerate these
        filter!(
            w -> any(j -> w[j] && !iszero(mags[j]), eachindex(mags)),
            nullspace,
        )
        patterns = [z]
        for w in nullspace
            append!(patterns, [p .⊻ w for p in patterns])
        end
    end
    x = MP.variables(model.model)
    for pattern in patterns
        values = zeros(T, length(x))
        for (j, var) in enumerate(vars)
            values[findfirst(isequal(var), x)] = pattern[j] ? -mags[j] : mags[j]
        end
        push!(
            solutions,
            PolyJuMP.Solution(values, model.model, model.feasibility_tolerance),
        )
    end
    return solutions
end
