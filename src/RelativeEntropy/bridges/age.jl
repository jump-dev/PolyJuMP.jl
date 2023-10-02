"""
    AGEBridge{T,F,G,H} <: MOI.Bridges.Constraint.AbstractBridge

The nonnegativity of `x ≥ 0` in
```
∑ ci x^αi ≥ -c0 x^α0
```
can be reformulated as
```
∑ ci exp(αi'y) ≥ -β exp(α0'y)
```
In any case, it is shown to be equivalent to
```
∃ ν ≥ 0 s.t. D(ν, exp(1)*c) ≤ β, ∑ νi αi = α0 ∑ νi [CP16, (3)]
```
where `N(ν, λ) = ∑ νj log(νj/λj)` is the relative entropy function.
The constant `exp(1)` can also be moved out of `D` into
```
∃ ν ≥ 0 s.t. D(ν, c) - ∑ νi ≤ β, ∑ νi αi = α0 ∑ νi [MCW21, (2)]
```
The relative entropy cone in MOI is `(u, v, w)` such that `D(w, v) ≥ u`.
Therefore, we can either encode `(β, exp(1)*c, ν)` [CP16, (3)] or
`(β + ∑ νi, c, ν)` [MCW21, (2)].
In this bridge, we use the second formulation.

!!! note
    A direct application of the Arithmetic-Geometric mean inequality gives
    ```
    ∃ ν ≥ 0 s.t. D(ν, exp(1)*c) ≤ -log(-β), ∑ νi αi = α0, ∑ νi = 1 [CP16, (4)]
    ```
    which is not jointly convex in (ν, c, β).
    The key to get the convex formulation [CP16, (3)] or [MCW21, (2)] instead is to
    use the convex conjugacy between the exponential and the negative entropy functions [CP16, (iv)].

[CP16] Chandrasekaran, Venkat, and Parikshit Shah.
"Relative entropy relaxations for signomial optimization."
SIAM Journal on Optimization 26.2 (2016): 1147-1173.
[MCW21] Murray, Riley, Venkat Chandrasekaran, and Adam Wierman.
"Signomial and polynomial optimization via relative entropy and partial dualization."
Mathematical Programming Computation 13 (2021): 257-295.
https://arxiv.org/pdf/1907.00814.pdf


"""
struct AGEBridge{T,F,G,H} <: MOI.Bridges.Constraint.AbstractBridge
    k::Int
    equality_constraints::Vector{MOI.ConstraintIndex{F,MOI.EqualTo{T}}}
    relative_entropy_constraint::MOI.ConstraintIndex{G,MOI.RelativeEntropyCone}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{AGEBridge{T,F,G,H}},
    model,
    func::G,
    set::AGECone,
) where {T,F,G,H}
    m = size(set.α, 1)
    ν = MOI.add_variables(model, m - 1)
    sumν = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(one(T), ν), zero(T))
    ceq = map(1:size(set.α, 2)) do var
        f = zero(typeof(sumν))
        j = 0
        for i in 1:m
            if i == set.k
                MA.sub_mul!!(f, convert(T, set.α[i, var]), sumν)
            else
                j += 1
                MA.add_mul!!(f, convert(T, set.α[i, var]), ν[j])
            end
        end
        return MOI.add_constraint(model, f, MOI.EqualTo(zero(T)))
    end
    scalars = MOI.Utilities.eachscalar(func)
    f = MOI.Utilities.operate(
        vcat,
        T,
        scalars[set.k] + sumν,
        scalars[setdiff(1:m, set.k)],
        MOI.VectorOfVariables(ν),
    )
    relative_entropy_constraint =
        MOI.add_constraint(model, f, MOI.RelativeEntropyCone(2m - 1))
    return AGEBridge{T,F,G,H}(set.k, ceq, relative_entropy_constraint)
end

function MOI.supports_constraint(
    ::Type{<:AGEBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:AGECone},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:AGEBridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:AGEBridge{T,F,G}},
) where {T,F,G}
    return [(F, MOI.EqualTo{T}), (G, MOI.RelativeEntropyCone)]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:AGEBridge{T}},
    H::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:AGECone},
) where {T}
    S = MOI.Utilities.scalar_type(H)
    F = MOI.Utilities.promote_operation(
        +,
        T,
        S,
        MOI.ScalarAffineFunction{Float64},
    )
    G = MOI.Utilities.promote_operation(vcat, T, S, F)
    return AGEBridge{T,F,G,H}
end
