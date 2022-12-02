module AlgebraicKKT

import MutableArithmetics
const MA = MutableArithmetics

import MathOptInterface
const MOI = MathOptInterface

import MultivariatePolynomials
const MP = MultivariatePolynomials

import SemialgebraicSets

import DynamicPolynomials

import PolyJuMP

const VarType = DynamicPolynomials.PolyVar{true}
const PolyType{T} = DynamicPolynomials.Polynomial{true,T}
const FuncType{T} = PolyJuMP.ScalarPolynomialFunction{T,PolyType}

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    # Model
    variables::Dict{MOI.VariableIndex,VarType}
    objective_sense::MOI.OptimizationSense
    objective_function::Union{Nothing,FuncType}
    set
    extrema::Vector{Vector{T}}
    objective_values::Vector{T}
    # Optimizer attributes
    algebraic_solver::Union{Nothing,SemialgebraicSets.AbstractAlgebraicSolver}
    feasibility_tolerance::T
    optimality_tolerance::T
    function Optimizer{T}() where {T}
        return new{T}(
            Dict{MOI.VariableIndex,VarType}(),
            MOI.FEASIBILITY_SENSE,
            nothing,
            SemialgebraicSets.FullSpace,
            Vector{T}[],
            T[],
            optimizer,
            Base.rtoldefault(T),
            Base.rtoldefault(T),
        )
    end
end

function MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, solver::SemialgebraicSets.AbstractAlgebraicSolver)
    if attr.name != "algebraic_solver"
        throw(MOI.UnsupportedAttribute(attr))
    end
    model.algebraic_solver = solver
end

function MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, tol)
    if attr.name == "feasibility_tolerance"
        model.feasibility_tolerance = tol
    elseif attr.name != "optimality_tolerance"
        model.optimality_tolerance = tol
    else
        throw(MOI.UnsupportedAttribute(attr))
    end
end

function MOI.is_empty(model::Optimizer)
    return isempty(model.variables) && model.objective_sense == MOI.FEASIBILITY_SENSE && isnothing(model.objective_function) && model.set isa SemialgebraicSets.FullSpace
end

function MOI.empty!(model::Optimizer)
    empty!(model.variables)
    model.objective_sense = MOI.FEASIBILITY_SENSE
    model.objective_function = nothing
    model.set = SemialgebraicSets.FullSpace()
    empty!(model.extrema)
    empty!(model.objective_values)
    return
end

function MOI.add_variable(model::Optimizer)
    i = length(model.variables) + 1
    vi = MOI.VariableIndex(i)
    var = VarType("x[$i]")
    model.variables[vi] = var
    return vi
end

function _polynomial(variables, func::PolyJuMP.ScalarPolynomialFunction)
    new_variables = VarType[variables[vi] for vi in func.variables]
    return func.polynomial(MP.variables(func.polynomial) => new_variables)
end

function _set(poly::MP.AbstractPolynomial{T}, set::MOI.EqualTo{T}) where {T}
    return SemialgebraicSets.equality(poly, MOI.constant(set))
end

function _set(poly::MP.AbstractPolynomial{T}, set::MOI.GreaterThan{T}) where {T}
    return SemialgebraicSets.PolynomialInequality(poly - MOI.constant(set))
end

function _set(poly::MP.AbstractPolynomial{T}, set::MOI.LessThan{T}) where {T}
    return SemialgebraicSets.PolynomialInequality(MOI.constant(set) - poly)
end

function MOI.add_constraint(model::Optimizer{T}, func::PolyJuMP.ScalarPolynomialFunction{T}, set::MOI.AbstractScalarSet) where {T}
    model.set = model.set ∩ _set(_polynomial(model.variables, func), set)
    return
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    model.objective_sense = sense
    return
end

function MOI.set(model::Optimizer{T}, ::MOI.ObjectiveFunction, func::PolyJuMP.ScalarPolynomialFunction{T}) where {T}
    model.objective_function = _polynomial(model.variables, func)
    return
end

function _add_to_system(system, lagrangian, set::SemialgebraicSets.FullSpace)
    return lagrangian
end

function _add_to_system(system, lagrangian, set::SemialgebraicSets.AlgebraicSet)
    n = SemialgebraicSets.nequalities(set)
    if iszero(n)
        return
    end
    DynamicPolynomials.@polyvar λ[1:n]
    for i in eachindex(λ)
        p = SemialgebraicSets.equalities(set)[i]
        SemialgebraicSets.addequality!(system, p)
        lagrangian = MA.add_mul!(lagrangian, λ[i], p)
    end
    return lagrangian
end

function _add_to_system(system, lagrangian, set::SemialgebraicSets.BasicSemialgebraicSet)
    lagrangian = _add_to_system(system, lagrangian, set.V)
    DynamicPolynomials.@polyvar σ[1:SemialgebraicSets.ninequalities(set)]
    for i in eachindex(σ)
        p = SemialgebraicSets.inequalities(set)[i]
        SemialgebraicSets.addequality!(system, σ[i] * p)
        lagrangian = MA.add_mul!(lagrangian, σ[i]^2, p)
    end
    return lagrangian
end

function _violates_inequalities(set::Union{SemialgebraicSets.FullSpace,SemialgebraicSets.AlgebraicSet}, x, sol, tol)
    return false
end

function _violates_inequalities(set::SemialgebraicSets.BasicSemialgebraicSet, x, sol, tol)
    return any(p -> p(x => sol) < -tol, SemialgebraicSets.inequalities(set))
end

function MOI.optimize!(model::Optimizer)
    DynamicPolynomials.@polyvar σ[1:ninequalities(model.set)]
    T = _coefficient_type(model)
    system = SemialgebraicSets.AlgebraicSet{T,polynomialtype(VarType, T)}()
    if model.objective_sense == MOI.FEASIBILITY_SENSE
        lagrangian = MA.Zero()
    else
        lagrangian = MA.mutable_copy(model.objective_function)
    end
    lagrangian = _add_to_system(system, lagrangian, model.set)
    x = collect(values(model.variables))
    ∇x = MP.differentiate(lagrangian, x)
    for p in ∇x
        SemialgebraicSets.addequality!(system, p)
    end
    model.extrema = Vector{T}[sol for sol in system if !_violates_inequalities(model.set, x, sol, model.tol)]
    if model.sense != MOI.FEASIBILITY_SENSE
        model.objective_values = T[model.objective_function(x => sol) for sol in model.extrema]
        I = sortperm(model.objective_values, rev = model.sense == MOI.MAX_SENSE)
        model.extrema = model.extrema[I]
        model.objective_values = model.objective_values[I]
        i = findfirst(eachindex(model.objective_values)) do i
            abs(model.objective_values[i] - first(model.objective_values)) > tol
        end
        if !isnothing(i)
            model.extrema = model.extrema[1:(i - 1)]
            model.objective_values = model.objective_values[1:(i - 1)]
        end
    end
    return
end

end
