struct SAGEBridge{T,F,G} <: MOI.Bridges.Constraint.AbstractBridge
    ν::Matrix{MOI.VariableIndex}
    age_constraints::Vector{MOI.ConstraintIndex{MOI.VectorOfVariables,SignomialAGECone}}
    equality_constraints::Vector{MOI.ConstraintIndex{F,MOI.EqualTo{T}}}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SAGEBridge{T,F,G}},
    model,
    func::G,
    set::SignomialSAGECone,
) where {T,F,G}
    m = size(set.α, 1)
    ν = Matrix{MOI.VariableIndex}(undef, m, m)
    age_constraints =
        Vector{MOI.ConstraintIndex{MOI.VectorOfVariables,SignomialAGECone}}(undef, m)
    for k in 1:m
        ν[k, :], age_constraints[k] =
            MOI.add_constrained_variables(model, SignomialAGECone(set.α, k))
    end
    scalars = MOI.Utilities.eachscalar(func)
    n = size(set.α, 2)
    equality_constraints = map(1:n) do i
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
    ::Type{<:SignomialSAGECone},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:SAGEBridge})
    return Tuple{Type}[(SignomialAGECone,)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:SAGEBridge{T,F}},
) where {T,F}
    return Tuple{Type,Type}[(F, MOI.EqualTo{T})]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SAGEBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SignomialSAGECone},
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
