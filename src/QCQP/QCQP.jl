module QCQP

import MutableArithmetics as MA
import MultivariatePolynomials as MP
import DataStructures
import PolyJuMP

# TODO special case for squares ?
function decompose(monos::AbstractVector{M}) where {M<:MP.AbstractMonomial}
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


function decompose(monos::AbstractVector{M}, ::Any) where {M<:MP.AbstractMonomial}
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

include("MOI_wrapper.jl")

end
