"""
    SignomialBridge{T,S,P,F} <: MOI.Bridges.Constraint.AbstractBridge

We use the Signomial Representative `SR` equation of [MCW21].

[MCW20] Riley Murray, Venkat Chandrasekaran, Adam Wierman
"Newton Polytopes and Relative Entropy Optimization"
https://arxiv.org/abs/1810.01614
[MCW21] Murray, Riley, Venkat Chandrasekaran, and Adam Wierman.
"Signomial and polynomial optimization via relative entropy and partial dualization."
Mathematical Programming Computation 13 (2021): 257-295.
https://arxiv.org/pdf/1907.00814.pdf
"""
struct SignomialBridge{T,S,P,F} <: MOI.Bridges.Constraint.AbstractBridge
    constraint::MOI.ConstraintIndex{F,S}
end

_signomial(set::PolynomialSAGECone) = SignomialSAGECone(set.α)
_signomial(::Type{PolynomialSAGECone}) = SignomialSAGECone
_signomial(set::PolynomialAGECone) = SignomialAGECone(set.α, set.k)
_signomial(::Type{PolynomialAGECone}) = SignomialAGECone

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SignomialBridge{T,S,P,F}},
    model,
    func::F,
    set,
) where {T,S,P,F}
    g = MOI.Utilities.scalarize(func)
    for i in eachindex(g)
        if any(isodd, set.α[i, :])
            vi = MOI.add_variable(model)
            # vi ≤ -|g[i]|
            MOI.Utilities.normalize_and_add_constraint(model, one(T) * vi - g[i], MOI.LessThan(zero(T)))
            MOI.Utilities.normalize_and_add_constraint(model, one(T) * vi + g[i], MOI.LessThan(zero(T)))
            g[i] = vi
        end
    end
    constraint = MOI.add_constraint(model, MOI.Utilities.vectorize(g), _signomial(set))
    return SignomialBridge{T,S,P,F}(constraint)
end

function MOI.supports_constraint(
    ::Type{<:SignomialBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:Union{PolynomialSAGECone,PolynomialAGECone}},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:SignomialBridge})
    return Tuple{Type}[(MOI.Reals,)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:SignomialBridge{T,S,P,F}},
) where {T,S,P,F}
    return [(F, S)]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SignomialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    P::Type{<:Union{PolynomialSAGECone,PolynomialAGECone}},
) where {T}
    S = _signomial(P)
    return SignomialBridge{T,S,P,F}
end

function MOI.get(
    model::MOI.ModelLike,
    attr::DecompositionAttribute,
    bridge::SignomialBridge,
)
    return MOI.get(model, attr, bridge.constraint)
end
