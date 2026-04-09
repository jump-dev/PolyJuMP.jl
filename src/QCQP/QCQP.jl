module QCQP

import MutableArithmetics as MA
import MultivariatePolynomials as MP
import DataStructures
import PolyJuMP

abstract type AbstractQuadratizationMethod end
struct GreedyQuadratization <: AbstractQuadratizationMethod end
struct BasicQuadratization <: AbstractQuadratizationMethod end
struct KACQuadratization <: AbstractQuadratizationMethod end

# TODO special case for squares ?
function decompose(monos::AbstractVector{M}, ::GreedyQuadratization) where {M<:MP.AbstractMonomial}
	println("RUNNING GREEDY QUADRATIZATION")
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


function decompose(monos::AbstractVector{M}, ::BasicQuadratization) where {M<:MP.AbstractMonomial}
	println("RUNNING BASIC QUADRATIZATION")
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
		for (var, n) in zip(mono.vars, mono.z)
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

# Inspired by https://www.sciencedirect.com/science/article/pii/S0098135422002472, by T. Karia, C. S. Adjiman, B. Chachuat => "KAC"
function decompose(monos::AbstractVector{M}, ::KACQuadratization) where {M<:MP.AbstractMonomial}
	println("RUNNING KAC QUADRATIZATION")
    vars = MP.variables(monos)
    quad = DataStructures.OrderedDict{M,Union{Nothing,M}}(
        var => var for var in vars
    )
	for x in vars
		for y in vars
			quad[x*y] = x
		end
	end

	stack = DataStructures.Stack{M}(length(monos))
	for mono in monos
		push!(stack, mono)
	end

	while !isempty(stack)
		mono = pop!(stack)
		deg = MP.degree(mono)
		if deg <= 2 || haskey(quad, mono)
			continue
		end

		stop = false
		for a in keys(quad)
			if MP.divides(a, mono)
				b = MP.div_multiple(mono, a)
				# println("$a divides $mono and gives $b")
				if haskey(quad, b)
					quad[mono] = a
					stop = true
					break
				elseif MP.degree(a) > 1
					quad[mono] = a
					push!(stack, b)
					stop = true
					break
				end
			end
		end
		if stop
			continue
		end

		l, r = nothing, nothing
		if maximum(mono.z) == 1
			v = filter(x -> MP.divides(x, mono), mono.vars)
			l = prod(v[1:deg÷2])
			r = prod(v[deg÷2+1:deg])
		elseif any(mono.z .% 2)
			l = MP.div_multiple(mono, mono)
			r = copy(mono)
			for i in 1:length(mono.vars)
				if r.z[i] % 2
					l *= r.vars[i]
					r.z[i] -= 1
				end
			end
		else
			l = copy(mono)
			l.z ./= 2
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


include("MOI_wrapper.jl")

end
