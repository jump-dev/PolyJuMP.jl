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

include("MOI_wrapper.jl")

end
