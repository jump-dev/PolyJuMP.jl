struct SAGEBridge{T,F,G} <: MOI.Bridges.Constraint.AbstractBridge
    ν::Matrix{MOI.VariableIndex}
    age_constraints::Vector{
        MOI.ConstraintIndex{MOI.VectorOfVariables,Cone{Signomials{Int}}},
    }
    equality_constraints::Vector{MOI.ConstraintIndex{F,MOI.EqualTo{T}}}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SAGEBridge{T,F,G}},
    model,
    func::G,
    set::Cone{Signomials{Nothing}},
) where {T,F,G}
    m = size(set.α, 1)
    ν = Matrix{MOI.VariableIndex}(undef, m, m)
    A = Cone{Signomials{Int}}
    age_constraints =
        Vector{MOI.ConstraintIndex{MOI.VectorOfVariables,A}}(undef, m)
    for k in 1:m
        ν[k, :], age_constraints[k] =
            MOI.add_constrained_variables(model, Cone(Signomials(k), set.α))
    end
    scalars = MOI.Utilities.eachscalar(func)
    equality_constraints = map(1:m) do i
        f = MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.(one(T), ν[:, i]),
            zero(T),
        )
        return MOI.Utilities.normalize_and_add_constraint(
            model,
            MA.sub!!(f, scalars[i]),
            MOI.EqualTo(zero(T)),
        )
    end
    return SAGEBridge{T,F,G}(ν, age_constraints, equality_constraints)
end

function MOI.supports_constraint(
    ::Type{<:SAGEBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{Cone{Signomials{Nothing}}},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:SAGEBridge})
    return Tuple{Type}[(Cone{Signomials{Int}},)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:SAGEBridge{T,F}},
) where {T,F}
    return Tuple{Type,Type}[(F, MOI.EqualTo{T})]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SAGEBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{Cone{Signomials{Nothing}}},
) where {T}
    S = MOI.Utilities.scalar_type(G)
    F = MOI.Utilities.promote_operation(
        -,
        T,
        MOI.ScalarAffineFunction{Float64},
        S,
    )
    return SAGEBridge{T,F,G}
end

function MOI.get(
    model::MOI.ModelLike,
    attr::DecompositionAttribute,
    bridge::SAGEBridge,
)
    return filter!(
        !isempty,
        [
            filter(
                x -> !isapprox(x, zero(x), atol=attr.tol),
                MOI.get(model, MOI.VariablePrimal(attr.result_index), bridge.ν[k, :])
            )
            for k in axes(bridge.ν, 1)
        ]
    )
end
