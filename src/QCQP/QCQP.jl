module QCQP

import MutableArithmetics as MA
import MultivariatePolynomials as MP
import DynamicPolynomials as DP
import DataStructures
import MathOptInterface as MOI
import PolyJuMP

"""
    abstract type AbstractQuadratizationMethod end

Method choosing how each monomial of degree larger than 2 is written as the
product of two monomials of lower degree, see [`decompose`](@ref) and the
[`QuadratizationMethod`](@ref) attribute.
"""
abstract type AbstractQuadratizationMethod end

"""
    GreedyQuadratization()

Monomials are decomposed by repeatedly adding the most frequent missing factor.
"""
struct GreedyQuadratization <: AbstractQuadratizationMethod end

"""
    @enum Decomposition Greedy Balanced Basic

Third phase of [`CompositeQuadratization`](@ref):

* `Greedy`: substitutes the most frequent pair of variables first;
* `Balanced`: splits multilinear monomials into two halves, requires the smart decomposition;
* `Basic`: sequential decomposition, separating the first two variables from the rest.
"""
@enum Decomposition begin
    Greedy
    Balanced
    Basic
end

"""
    CompositeQuadratization(;
        decompo::Decomposition = Balanced,
        smart::Bool = true,
        variable_reuse::Bool = true,
    )

Multi-level quadratization: if `smart`, the smart decomposition first exploits
even exponents and leaves only multilinear monomials, which are then decomposed
by `decompo`. With `variable_reuse`, a monomial that is the product of two
monomials already known is decomposed as such.
"""
struct CompositeQuadratization <: AbstractQuadratizationMethod
    decomposition::Decomposition
    smart::Bool
    variable_reuse::Bool
    function CompositeQuadratization(;
        decompo::Decomposition = Balanced,
        smart::Bool = true,
        variable_reuse::Bool = true,
    )
        if decompo == Balanced && !smart
            throw(ArgumentError("The balanced decomposition requires `smart = true`."))
        end
        return new(decompo, smart, variable_reuse)
    end
end

"""
    OptimalQuadratization(optimizer::MOI.AbstractOptimizer)

Minimizes the number of auxiliary variables by solving a binary program with
`optimizer`. The candidate auxiliary variables are the products of pairs of
variables and products of these candidates, so the result is minimal among these candidates.
"""
struct OptimalQuadratization <: AbstractQuadratizationMethod
    optimizer::MOI.AbstractOptimizer
end


# TODO special case for squares ?
function decompose(monos::AbstractVector{M}, ::GreedyQuadratization; verbose::Bool=true) where {M<:MP.AbstractMonomial}
    if verbose
        println("Running greedy quadratization scheme...")
    end
    vars = MP.variables(monos)
    quad = DataStructures.OrderedDict{M,Union{Nothing,M}}(
        var => var for var in vars
    )
    for mono in monos
        if MP.degree(mono) > 1
            quad[mono] = nothing
        end
    end
    candidates =
        DataStructures.PriorityQueue{eltype(monos),Int}(Base.Order.Reverse)
    while any(mono -> isnothing(quad[mono]), keys(quad))
        empty!(candidates)
        for mono in keys(quad)
            if !isnothing(quad[mono])
                continue
            end
            for a in keys(quad)
                if a != mono && MP.divides(a, mono)
                    b = MP.div_multiple(mono, a)
                    if b in keys(quad)
                        quad[mono] = a
                    else
                        candidates[b] = get(candidates, b, 0) + 1
                    end
                end
            end
        end
        if !isempty(candidates)
            next, _ = first(candidates)
            quad[next] = nothing
        end
    end
    return quad
end


"""
    kac_phase(monos; variable_reuse::Bool, multilinear::Bool)

Decomposition inspired by Karia, Adjiman and Chachuat (2022), shared by the
smart decomposition (`multilinear = false`) and the
balanced decomposition (`multilinear = true`): monomials with an
odd exponent are split into their odd part and their even part, monomials
with even exponents into the square of their half and, if `multilinear`,
multilinear monomials into two halves; otherwise multilinear monomials are
left undecomposed (`nothing`) for the next phase.
"""
function kac_phase(monos::AbstractVector{M}; variable_reuse::Bool, multilinear::Bool, verbose::Bool=true) where {M<:MP.AbstractMonomial}
    # This function exists because smart decomposition into balanced decomposition
    # used to be one method, so the code evolved into having this function, and when
    # the method was split this function was kept but parametrized to avoid breaking anything,
    # it probably could do with a refactor at some point
    vars = MP.variables(monos)
    quad = DataStructures.OrderedDict{M,Union{Nothing,M}}()

    stack = DataStructures.Stack{M}()
    for mono in monos
        if MP.degree(mono) == 0
            continue
        end
        quad[mono] = nothing
        push!(stack, mono)
    end

    while !isempty(stack)
        mono = pop!(stack)
        deg = MP.degree(mono)
        if deg <= 2 || get(quad, mono, nothing) !== nothing
            continue
        end

        if variable_reuse
            stop = false
            for a in keys(quad)
                if a != mono && MP.divides(a, mono)
                    b = MP.div_multiple(mono, a)
                    if haskey(quad, b)
                        quad[mono] = a
                        stop = true
                        break
                    end
                end
            end
            if stop
                continue
            end
        end

        l, r = nothing, nothing
        mono_vars = MP.variables(mono)
        exps = MP.exponents(mono)
        if maximum(exps) == 1 && multilinear
            v = filter(x -> MP.degree(mono, x) == 1, mono_vars)
            l = prod(v[1:deg÷2])
            r = prod(v[deg÷2+1:deg])
        elseif maximum(exps) == 1 && !multilinear
            quad[mono] = nothing
            continue
        elseif any((exps .% 2) .== 1)
            l = MP.div_multiple(mono, mono)
            r = copy(mono)
            for (v, e) in MP.powers(mono)
                if e % 2 == 1
                    l *= v
                    r = MP.div_multiple(r, v)
                end
            end
        else
            l = MP.div_multiple(mono, mono)
            for (v, e) in MP.powers(mono)
                l *= v^(e÷2)
            end
        end

        if l !== nothing
            quad[mono] = l
            push!(stack, l)
        end
        if r !== nothing
            push!(stack, r)
        end
    end

    return quad
end


function smart_decompo_phase(monos::AbstractVector{M}; variable_reuse::Bool, verbose::Bool=true) where {M<:MP.AbstractMonomial}
    return kac_phase(monos, variable_reuse=variable_reuse, multilinear=false, verbose=verbose)
end


function cut_in_half_phase(monos::AbstractVector{M}; variable_reuse::Bool, verbose::Bool=true) where {M<:MP.AbstractMonomial}
    return kac_phase(monos, variable_reuse=variable_reuse, multilinear=true, verbose=verbose)
end


# `greedy_phase` creates new polynomial variables for the selected pairs. With
# TypedPolynomials the variables are part of the monomial type, so new
# variables cannot be created at runtime (this used to throw a
# `MethodError`). The phase is therefore run on a DynamicPolynomials copy of
# the monomials, with the variables created in the same order so that the
# monomial order (and thus the tie-breaking) is the same, and the result is
# mapped back. DynamicPolynomials monomials, used by `QCQP.Optimizer`, take
# the direct path below, unchanged.
function greedy_phase(monos::AbstractVector{M}; verbose::Bool=true) where {M<:MP.AbstractMonomial}
    vars = MP.variables(monos)
    if isempty(vars)
        return DataStructures.OrderedDict{M,Union{Nothing,M}}()
    end
    DP.@polyvar dvars[1:length(vars)]
    to_dynamic(m) = prod(dvars[i]^MP.degree(m, vars[i]) for i in eachindex(vars); init = MP.constant_monomial(dvars[1]))
    from_dynamic(m) = convert(M, prod(vars[i]^MP.degree(m, dvars[i]) for i in eachindex(vars); init = MP.constant_monomial(vars[1])))
    dquad = greedy_phase(map(to_dynamic, monos); verbose = verbose)
    quad = DataStructures.OrderedDict{M,Union{Nothing,M}}()
    for (k, v) in dquad
        quad[from_dynamic(k)] = v === nothing ? nothing : from_dynamic(v)
    end
    return quad
end

function greedy_phase(monos::AbstractVector{M}; verbose::Bool=true) where {M<:DP.Monomial}
    vars = MP.variables(monos)
    quad = DataStructures.OrderedDict{M,Union{Nothing,M}}(
        var => var for var in vars
    )
    for x in vars
        for y in vars
            quad[x*y] = x
        end
    end

    var_to_mono = Dict{M, M}(var => var for var in vars)
    mono_to_var = Dict{M, M}()

    occurrence = DataStructures.PriorityQueue{M, Tuple{Int, M}}(Base.Order.Reverse)

    aux_var_index = 1
    monos_list = M[mono for mono in monos]
    while any(x -> MP.degree(x) > 2, monos_list)
        empty!(occurrence)
        for mono in monos_list
            if MP.degree(mono) <= 2
                continue
            end
            vs = MP.variables(mono)
            exps = MP.exponents(mono)
            nonzeroidx = [i for i in eachindex(exps) if exps[i] > 0]
            for k in eachindex(nonzeroidx)
                i = nonzeroidx[k]
                for j in nonzeroidx[k:end]
                    if i == j && exps[i] <= 1
                        continue
                    end
                    q = vs[i]*vs[j]
                    curr, _ = get(occurrence, q, (0, q))
                    occ = i == j ? exps[i]÷2 : min(exps[i], exps[j])
                    occurrence[q] = (curr + occ, q)
                end
            end
        end

        for (i, mono) in enumerate(monos_list)
            if MP.degree(mono) <= 2
                continue
            end
            for a in keys(occurrence)
                if MP.divides(a, mono)
                    deg = 0
                    b = mono
                    while MP.divides(a, b)
                        deg += 1
                        b = MP.div_multiple(b, a)
                    end
                    if haskey(mono_to_var, a)
                        v = mono_to_var[a]
                    else
                        DP.@polyvar v[aux_var_index]
                        aux_var_index += 1
                        v = v[1]
                        var_to_mono[v] = a
                        mono_to_var[a] = v
                    end
                    monos_list[i] = b*v^deg
                    break
                end
            end
        end
    end

    _evaluate(x) = begin
        mono = var_to_mono[x]
        vars = MP.variables(mono)
        if MP.degree(mono) == 1
            return x
        else
            nonzero_vars = filter(v -> MP.degree(mono, v) >= 1, vars)
            if length(nonzero_vars) == 1
                nonzero_vars = [nonzero_vars[1], nonzero_vars[1]]
            end
            return prod(_evaluate.(nonzero_vars), init=one(M))
        end
    end
    
    stack = DataStructures.Stack{M}()
    for mono in monos_list
        if MP.degree(mono) > 1
            push!(stack, mono)
        end
    end

    while !isempty(stack)
        bilin = pop!(stack)
        vars = filter(v -> MP.degree(bilin, v) >= 1, MP.variables(bilin))
        a, b = length(vars) == 2 ? (vars[1], vars[2]) : (vars[1], vars[1])
        a_mono, b_mono = _evaluate(a), _evaluate(b)
        mono = a_mono * b_mono

        if haskey(quad, mono) && !isnothing(quad[mono])
            continue
        end

        quad[mono] = a_mono

        if !haskey(quad, a_mono) && MP.degree(a_mono) > 2
            quad[a_mono] = nothing
            push!(stack, var_to_mono[a])
        end
        if !haskey(quad, b_mono) && MP.degree(b_mono) > 2
            quad[b_mono] = nothing
            push!(stack, var_to_mono[b])
        end
    end

    return quad
end


function fallback_phase(monos::AbstractVector{M}; verbose::Bool=true) where {M<:MP.AbstractMonomial}
    vars = MP.variables(monos)
    quad = DataStructures.OrderedDict{M,Union{Nothing,M}}(
        var => var for var in vars
    )
    for x in vars
        for y in vars
            quad[x*y] = x
        end
    end

    for mono in monos
        deg = MP.degree(mono)
        if deg <= 2
            continue
        end

        terms = M[]
        for (var, n) in MP.powers(mono)
            for _ in 1:n
                push!(terms, var)
            end
        end

        a = (deg-1)÷2
        epp = prod(terms[2*a+1:deg])
        quad[mono] = epp
        
        terms2 = terms[1:2*a]
        ep = prod(terms2)
        while MP.degree(ep) >= 2
            if haskey(quad, ep)
                break
            end
            div = prod(terms2[1:2])
            quad[ep] = div
            terms2 = terms2[3:end]
            ep = prod(terms2)
        end
    end

    return quad
end


function decompose(monos::AbstractVector{M}, method::OptimalQuadratization; verbose::Bool=true) where {M<:MP.AbstractMonomial}
    if verbose
        println("Running minimum auxiliary variables quadratization scheme...")
    end

    vars = MP.variables(monos)
    quad = DataStructures.OrderedDict{M,Union{Nothing,M}}(
        var => var for var in vars
    )
    deg = maximum(MP.degree.(monos))

    var_to_mono = DataStructures.OrderedDict{MOI.VariableIndex, M}()
    var_divisor = DataStructures.OrderedDict{MOI.VariableIndex, Union{M,MOI.VariableIndex}}()

    optimizer = method.optimizer
    monos = filter(x -> MP.degree(x) >= 3, monos)
    tmp = Set{M}()
    # find all pairs of variables
    for mono in monos
        vs = MP.variables(mono)
        exps = MP.exponents(mono)
        nonzeroidx = [i for i in eachindex(exps) if exps[i] > 0]
        for k in eachindex(nonzeroidx)
            i = nonzeroidx[k]
            for j in nonzeroidx[k:end]
                bilin = vs[i]*vs[j]
                if (i == j && exps[i] <= 1) || bilin in tmp
                    continue
                end
                new_var, _ = MOI.add_constrained_variable(optimizer, MOI.ZeroOne())
                var_to_mono[new_var] = bilin
                var_divisor[new_var] = MP.div_multiple(bilin, vs[j])
                push!(tmp, bilin)
            end
        end
    end
    empty!(tmp)

    if verbose
        println("found all pairs...")
    end

    explored_pairs = Set{Tuple{MOI.VariableIndex,Union{MOI.VariableIndex,eltype(vars)}}}()
    curr_deg = 2
    while curr_deg < deg-1
        aux_vars = collect(keys(var_to_mono))
        for (i, va) in enumerate(aux_vars)
            a = var_to_mono[va]
            for vb in aux_vars[i:end]
                b = var_to_mono[vb]
                if (va, vb) in explored_pairs || (vb, va) in explored_pairs
                    continue
                end
                push!(explored_pairs, (va, vb))
                stop = true
                for mono in monos
                    if MP.divides(a*b, mono)
                        stop = false
                        break
                    end
                end
                if stop
                    continue
                end

                if MP.degree(a*b) >= deg
                    continue
                end
                new_var, _ = MOI.add_constrained_variable(optimizer, MOI.ZeroOne())

                x = MOI.ScalarAffineTerm(1.0, new_var)
                y = MOI.ScalarAffineTerm(-1.0, va)
                MOI.add_constraint(optimizer, MOI.ScalarAffineFunction([x, y], 0.0), MOI.LessThan(0.0))

                x = MOI.ScalarAffineTerm(1.0, new_var)
                y = MOI.ScalarAffineTerm(-1.0, vb)
                MOI.add_constraint(optimizer, MOI.ScalarAffineFunction([x, y], 0.0), MOI.LessThan(0.0))

                var_to_mono[new_var] = a*b
                var_divisor[new_var] = va
            end
            for v in vars
                if (va, v) in explored_pairs
                    continue
                end
                push!(explored_pairs, (va, v))
                stop = true
                for mono in monos
                    if MP.divides(a*v, mono)
                        stop = false
                        break
                    end
                end
                if stop
                    continue
                end

                new_var, _ = MOI.add_constrained_variable(optimizer, MOI.ZeroOne())

                x = MOI.ScalarAffineTerm(1.0, new_var)
                y = MOI.ScalarAffineTerm(-1.0, va)
                MOI.add_constraint(optimizer, MOI.ScalarAffineFunction([x, y], 0.0), MOI.LessThan(0.0))

                var_to_mono[new_var] = a*v
                var_divisor[new_var] = va
            end
        end
        if verbose
            println("iteration deg = $curr_deg out of $deg done...")
        end
        curr_deg *= 2
    end

    if verbose
        println("done with variable creation...")
    end

    mono_to_vars = Dict{M, Vector{MOI.VariableIndex}}()
    for (v, m) in var_to_mono
        if !haskey(mono_to_vars, m)
            mono_to_vars[m] = MOI.VariableIndex[]
        end
        push!(mono_to_vars[m], v)
    end

    for mono in monos
        quad_terms = MOI.ScalarQuadraticTerm{Float64}[]
        affine_terms = MOI.ScalarAffineTerm{Float64}[]
        for (a, vs) in mono_to_vars
            if MP.divides(a, mono)
                b = MP.div_multiple(mono, a)
                if b == a || MP.degree(b) == 1
                    for v in vs
                        push!(affine_terms, MOI.ScalarAffineTerm(1.0, v))
                    end
                elseif MP.degree(b) > 1
                    for va in vs
                        for vb in mono_to_vars[b]
                            push!(quad_terms, MOI.ScalarQuadraticTerm(1.0, va, vb))
                        end
                    end
                end
            end
        end
        MOI.add_constraint(optimizer, MOI.ScalarQuadraticFunction(quad_terms, affine_terms, 0.0), MOI.GreaterThan(1.0))
    end

    if verbose
        println("done with constraint creation...")
    end

    affine_terms = MOI.ScalarAffineTerm{Float64}[]
    for var in keys(var_to_mono)
        push!(affine_terms, MOI.ScalarAffineTerm(1.0, var))
    end
    obj = MOI.ScalarAffineFunction(affine_terms, 0.0)
    MOI.set(optimizer, MOI.ObjectiveFunction{typeof(obj)}(), obj)
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    MOI.optimize!(optimizer)
    status = MOI.get(optimizer, MOI.TerminationStatus())
    if status != MOI.OPTIMAL
        println("OptimalQuadratization: the binary program selecting the auxiliary variables terminated with status $status.")
    end

    count_new_vars = 0
    active_vars = MOI.VariableIndex[]
    for (var, mono) in var_to_mono
        active = MOI.get(optimizer, MOI.VariablePrimal(), var)
        if active ≈ 1
            if verbose
                println("$mono is active")
            end
            if haskey(quad, mono) && verbose
                # This should never happen
                println("[WARNING] quad already has $mono, which means that $mono was activated twice...")
            end
            if var_divisor[var] isa M
                quad[mono] = var_divisor[var]
            else
                quad[mono] = var_to_mono[var_divisor[var]]
            end
            count_new_vars += 1
            push!(active_vars, var)
        end
    end

    if verbose
        println("should generate $count_new_vars new variables")
    end
    
    mono_set = Set{M}(monos)
    for (i, va) in enumerate(active_vars)
        a = var_to_mono[va]
        for vb in active_vars[i:end]
            b = var_to_mono[vb]
            mono = a*b
            if mono in mono_set && !haskey(quad, mono)
                quad[mono] = a
            end
        end
    end
    for mono in monos
        if haskey(quad, mono)
            continue
        end
        for (var, exp) in MP.powers(mono)
            if exp > 0 && haskey(quad, MP.div_multiple(mono, var))
                quad[mono] = MP.div_multiple(mono, var)
                break
            end
        end
    end

    for x in vars
        for y in vars
            quad[x*y] = x
        end
    end

    for mono in keys(quad)
        if quad[mono] === nothing
            error("OptimalQuadratization: no decomposition was selected for $mono.")
        end
    end

    MOI.empty!(optimizer)

    return quad
end


function decompose(monos::AbstractVector{M}, method::CompositeQuadratization; verbose::Bool=true) where {M<:MP.AbstractMonomial}
    if verbose
        println("Running composite quadratization scheme...")
        if method.smart
            print("Smart ")
        end
        print(lowercase(string(method.decomposition)))
        if method.variable_reuse
            print(" with variable reuse")
        end
        print("\n")
    end

    vars = MP.variables(monos)
    quad = DataStructures.OrderedDict{M,Union{Nothing,M}}()

    if method.smart
        quad = smart_decompo_phase(monos, variable_reuse=method.variable_reuse)
        monos = filter(x -> isnothing(quad[x]), collect(keys(quad)))
    end

    if method.decomposition == Greedy
        q = greedy_phase(monos, verbose=verbose)
    elseif method.decomposition == Balanced  # implies `method.smart`
        q = cut_in_half_phase(monos, variable_reuse=method.variable_reuse, verbose=verbose)
    else
        q = fallback_phase(monos, verbose=verbose)
    end

    for (k, v) in q
        if isnothing(get(quad, k, nothing))
            quad[k] = v
        end
    end

    for x in vars
        quad[x] = x
        for y in vars
            quad[x*y] = x
        end
    end

    return quad
end


include("MOI_wrapper.jl")

end
