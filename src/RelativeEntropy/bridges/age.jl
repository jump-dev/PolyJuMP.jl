struct AGEBridge{T,F,G,H} <: MOI.Bridges.Constraint.AbstractBridge
    k::Int
    equality_constraints::Vector{MOI.ConstraintIndex{F,MOI.EqualTo{T}}}
    relative_entropy_constraint::MOI.ConstraintIndex{G,MOI.RelativeEntropyCone}
end # See https://arxiv.org/pdf/1907.00814.pdf

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
        f = -sumν
        j = 0
        for i in 1:m
            if i != set.k
                j += 1
                MA.add_mul!!(f, convert(T, set.α[i, var]), ν[j])
            end
        end
        return MOI.Utilities.normalize_and_add_constraint(
            model,
            f,
            MOI.EqualTo(zero(T)),
        )
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
