"""
    struct ScalarPolynomialFunction{T,P<:MP.AbstractPolynomial{T}} <: MOI.AbstractScalarFunction
        polynomial::P
        variables::Vector{MOI.VariableIndex}
    end

Defines the polynomial function of the variables `variables` where the variable
`variables(p)[i]` corresponds to `variables[i]`.
"""
struct ScalarPolynomialFunction{T,P<:MP.AbstractPolynomial{T}} <:
       MOI.AbstractScalarFunction
    polynomial::P
    variables::Vector{MOI.VariableIndex}
end

function MOI.constant(func::ScalarPolynomialFunction)
    return MP.coefficient(
        func.polynomial,
        MP.constant_monomial(func.polynomial),
    )
end

function _polynomial_variable(::Type{P}, vi::MOI.VariableIndex) where {P}
    return MP.similar_variable(P, Symbol("x[$(vi.value)]"))
end

function Base.convert(
    ::Type{ScalarPolynomialFunction{T,P}},
    vi::MOI.VariableIndex,
) where {T,P}
    x = _polynomial_variable(P, vi)
    return ScalarPolynomialFunction{T,P}(MP.polynomial(x, T), [vi])
end

function _polynomial_variables!(::Type{P}, variables) where {P}
    sort!(variables, by = v -> v.value)
    unique!(variables)
    # FIXME It is an issue for TypedPolynomials since `"x[1]" < "x[2]"` but `"x[1]" > "x[10]"`
    x = _polynomial_variable.(P, variables)
    if !issorted(x, rev = true)
        error("`$P` unsupported, use DynamicPolynomials instead")
    end
    d = Dict(variables[i] => x[i] for i in eachindex(variables))
    return x, d
end

function _polynomial_with_variables(
    ::Type{P},
    func::MOI.ScalarAffineFunction,
    d,
) where {P}
    terms = [MP.term(t.coefficient, d[t.variable]) for t in func.terms]
    push!(terms, MOI.constant(func))
    return MP.polynomial(terms)
end

function Base.convert(
    ::Type{ScalarPolynomialFunction{T,P}},
    func::MOI.ScalarAffineFunction{T},
) where {T,P}
    variables = [t.variable for t in func.terms]
    _, d = _polynomial_variables!(P, variables)
    poly = _polynomial_with_variables(P, func, d)
    return ScalarPolynomialFunction{T,P}(poly, variables)
end

function _to_polynomial!(
    d::Dict{K,V},
    ::Type{T},
    func::MOI.ScalarQuadraticFunction{T},
) where {K,V,T}
    terms = MP.term_type(V, T)[MOI.constant(func)]
    for t in func.affine_terms
        push!(terms, MP.term(t.coefficient, _to_polynomial!(d, T, t.variable)))
    end
    for t in func.quadratic_terms
        coef = t.variable_1 == t.variable_2 ? t.coefficient / 2 : t.coefficient
        push!(
            terms,
            MP.term(
                coef,
                _to_polynomial!(d, T, t.variable_1) *
                _to_polynomial!(d, T, t.variable_2),
            ),
        )
    end
    return MP.polynomial(terms)
end

function Base.convert(
    ::Type{ScalarPolynomialFunction{T,P}},
    func::MOI.ScalarQuadraticFunction{T},
) where {T,P}
    linear_variables = [t.variable for t in func.affine_terms]
    quad_variables_1 = [t.variable_1 for t in func.quadratic_terms]
    quad_variables_2 = [t.variable_2 for t in func.quadratic_terms]
    variables = [linear_variables; quad_variables_1; quad_variables_2]
    _, d = _polynomial_variables!(P, variables)
    poly = _to_polynomial!(d, T, func)
    return ScalarPolynomialFunction{T,P}(poly, variables)
end

function Base.convert(
    ::Type{ScalarPolynomialFunction{T,P}},
    func::MOI.ScalarNonlinearFunction,
) where {T,P}
    return _to_polynomial(func, T)
end

function Base.copy(func::ScalarPolynomialFunction)
    return ScalarPolynomialFunction(
        MA.copy_if_mutable(func.polynomial),
        copy(func.variables),
    )
end

function MOI.Utilities.canonicalize!(::ScalarPolynomialFunction) end

function _variables(aff::MOI.ScalarAffineFunction)
    return MOI.VariableIndex[t.variable for t in aff.terms]
end

function MOI.Utilities.substitute_variables(
    variable_map::Function,
    func::ScalarPolynomialFunction{T,P},
) where {T,P}
    new_aff =
        MOI.ScalarAffineFunction{T}[variable_map(var) for var in func.variables]
    variables = collect(Iterators.flatten(_variables(aff) for aff in new_aff))
    _, d = _polynomial_variables!(P, variables)
    new_polys = [_polynomial_with_variables(P, aff, d) for aff in new_aff]
    new_poly = func.polynomial(MP.variables(func.polynomial) => new_polys)
    return ScalarPolynomialFunction{T,typeof(new_poly)}(new_poly, variables)
end

function MOI.Utilities.is_coefficient_type(
    ::Type{<:ScalarPolynomialFunction{T}},
    ::Type{S},
) where {S,T}
    return S === T
end

# Placeholder for `promote_operation`
struct VectorPolynomialFunction{T,P<:MP.AbstractPolynomial{T}} <:
       MOI.AbstractVectorFunction end

function MOI.Utilities.scalar_type(
    ::Type{VectorPolynomialFunction{T,P}},
) where {T,P}
    return PolyJuMP.ScalarPolynomialFunction{T,P}
end

function MOI.Utilities.is_coefficient_type(
    ::Type{<:VectorPolynomialFunction{T}},
    ::Type{T},
) where {T}
    return true
end

function MOI.Utilities.is_coefficient_type(
    ::Type{<:VectorPolynomialFunction},
    ::Type,
)
    return false
end

function MOI.Utilities.promote_operation(
    ::typeof(-),
    ::Type{T},
    F::Type{
        <:Union{ScalarPolynomialFunction{T,P},VectorPolynomialFunction{T,P}},
    },
) where {T,P}
    return F
end

function MOI.Utilities.promote_operation(
    ::typeof(-),
    ::Type{T},
    F::Type{ScalarPolynomialFunction{T,P}},
    ::Type{<:Union{T,MOI.Utilities.ScalarLike{T}}},
) where {T,P}
    return F
end

function MOI.Utilities.promote_operation(
    ::typeof(-),
    ::Type{T},
    F::Type{VectorPolynomialFunction{T,P}},
    ::Type{<:Union{AbstractVector{T},MOI.Utilities.VectorLike{T}}},
) where {T,P}
    return F
end

function MOI.Utilities.promote_operation(
    ::typeof(vcat),
    ::Type{T},
    ::Type{ScalarPolynomialFunction{T,P}},
) where {T,P}
    return VectorPolynomialFunction{T,P}
end

function MOI.Utilities.operate(
    op::Union{typeof(+),typeof(-)},
    ::Type{T},
    p::ScalarPolynomialFunction{T,P},
    f::Union{T,MOI.AbstractScalarFunction},
) where {T,P}
    d = Dict(
        vi => v for (vi, v) in zip(p.variables, MP.variables(p.polynomial))
    )
    poly = _to_polynomial!(d, T, f)
    return _scalar_polynomial(d, T, op(p.polynomial, poly))
end
