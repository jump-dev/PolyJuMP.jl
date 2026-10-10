"""
    SignomialsBridge{T,S,P,F} <: MOI.Bridges.Constraint.AbstractBridge

We use the Signomials Representative `SR` equation of [MCW21].

[MCW20] Riley Murray, Venkat Chandrasekaran, Adam Wierman
"Newton Polytopes and Relative Entropy Optimization"
https://arxiv.org/abs/1810.01614
[MCW21] Murray, Riley, Venkat Chandrasekaran, and Adam Wierman.
"Signomials and polynomial optimization via relative entropy and partial dualization."
Mathematical Programming Computation 13 (2021): 257-295.
https://arxiv.org/pdf/1907.00814.pdf
"""
struct SignomialsBridge{T,S,P,F,G} <: MOI.Bridges.Constraint.AbstractBridge
    # Indices of the rows `i` of `set.α` that have an odd entry
    odd::Vector{Int}
    # For each odd row `i`, the constraint `vi - g[i] ≤ 0`
    lower::Vector{MOI.ConstraintIndex{G,MOI.LessThan{T}}}
    # For each odd row `i`, the constraint `vi + g[i] ≤ 0`
    upper::Vector{MOI.ConstraintIndex{G,MOI.LessThan{T}}}
    constraint::MOI.ConstraintIndex{F,S}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SignomialsBridge{T,S,P,F,G}},
    model,
    func::F,
    set,
) where {T,S,P,F,G}
    g = MOI.Utilities.scalarize(func)
    odd = Int[]
    lower = MOI.ConstraintIndex{G,MOI.LessThan{T}}[]
    upper = MOI.ConstraintIndex{G,MOI.LessThan{T}}[]
    for i in eachindex(g)
        if any(isodd, set.α[i, :])
            vi = MOI.add_variable(model)
            push!(odd, i)
            # vi ≤ -|g[i]|
            push!(
                lower,
                MOI.Utilities.normalize_and_add_constraint(
                    model,
                    one(T) * vi - g[i],
                    MOI.LessThan(zero(T)),
                ),
            )
            push!(
                upper,
                MOI.Utilities.normalize_and_add_constraint(
                    model,
                    one(T) * vi + g[i],
                    MOI.LessThan(zero(T)),
                ),
            )
            g[i] = vi
        end
    end
    constraint = MOI.add_constraint(
        model,
        MOI.Utilities.vectorize(g),
        Cone(Signomials(set.cone.monomial), set.α),
    )
    return SignomialsBridge{T,S,P,F,G}(odd, lower, upper, constraint)
end

function MOI.supports_constraint(
    ::Type{<:SignomialsBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{Cone{Polynomials{M}}},
) where {T,M}
    return true
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:SignomialsBridge},
)
    return Tuple{Type}[(MOI.Reals,)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:SignomialsBridge{T,S,P,F,G}},
) where {T,S,P,F,G}
    return [(F, S), (G, MOI.LessThan{T})]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SignomialsBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    P::Type{Cone{Polynomials{M}}},
) where {T,M}
    G = MOI.Utilities.promote_operation(
        -,
        T,
        MOI.ScalarAffineFunction{T},
        MOI.Utilities.scalar_type(F),
    )
    return SignomialsBridge{T,Cone{Signomials{M}},P,F,G}
end

function MOI.get(
    model::MOI.ModelLike,
    attr::DecompositionAttribute,
    bridge::SignomialsBridge,
)
    return MOI.get(model, attr, bridge.constraint)
end

# The dual of the polynomial SAGE constraint is the adjoint of the linear map
# used in the reformulation, applied to the duals of the constraints created
# by the bridge. For a row `i` with even exponents, `g[i]` only appears in row
# `i` of the signomial constraint so the dual is the corresponding entry `w[i]`
# of its dual `w`. For an odd row, `g[i]` appears with coefficient `-1` in
# `lower[i]` and `+1` in `upper[i]` so the dual is the difference of the duals
# of these two constraints. The result `v` satisfies `|v[i]| ≤ w[i]` for odd
# rows, which matches the characterization of the dual of the polynomial SAGE
# cone in terms of the dual of the signomial SAGE cone of [MCW20].
function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::SignomialsBridge,
)
    v = MOI.get(model, attr, bridge.constraint)
    for (j, i) in enumerate(bridge.odd)
        v[i] =
            MOI.get(model, attr, bridge.upper[j]) -
            MOI.get(model, attr, bridge.lower[j])
    end
    return v
end
