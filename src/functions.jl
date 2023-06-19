"""
    struct ScalarPolynomialFunction{P<:MP.AbstractPolynomialLike} <: MOI.AbstractScalarFunction
        p::P
        variables::Vector{MOI.VariableIndex}
    end

Defines the polynomial function of the variables `variables` where the variable
`variables(p)[i]` corresponds to `variables[i]`.
"""
struct ScalarPolynomialFunction{T,P<:AbstractPolynomial{T}} <:
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
    return ScalarPolynomialFunction{T,P}(polynomial(x, T), [vi])
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
    x, d = _polynomial_variables!(P, variables)
    poly = _polynomial_with_variables(P, func, d)
    return ScalarPolynomialFunction{T,P}(poly, variables)
end

function Base.convert(
    ::Type{ScalarPolynomialFunction{T,P}},
    func::MOI.ScalarQuadraticFunction{T},
) where {T,P}
    linear_variables = [t.variable for t in func.affine_terms]
    quad_variables_1 = [t.variable_1 for t in func.quadratic_terms]
    quad_variables_2 = [t.variable_2 for t in func.quadratic_terms]
    variables = [linear_variables; quad_variables_1; quad_variables_2]
    x, d = _polynomial_variables!(P, variables)
    terms = MP.term_type(P)[MOI.constant(func)]
    for t in func.affine_terms
        push!(terms, MP.term(t.coefficient, d[t.variable]))
    end
    for t in func.quadratic_terms
        coef = t.variable_1 == t.variable_2 ? t.coefficient / 2 : t.coefficient
        push!(terms, MP.term(coef, d[t.variable_1] * d[t.variable_2]))
    end
    return ScalarPolynomialFunction{T,P}(polynomial(terms), variables)
end

function Base.convert(
    ::Type{FuncType{T}},
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
    x, d = _polynomial_variables!(P, variables)
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

function MOI.Utilities.promote_operation(
    ::typeof(-),
    ::Type{T},
    F::Type{ScalarPolynomialFunction{T,P}},
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

# FIXME
function MOI.Utilities.promote_operation(
    ::typeof(vcat),
    ::Type{T},
    ::Type{ScalarPolynomialFunction{T,P}},
) where {T,P}
    return MOI.VectorQuadraticFunction{T}
end
