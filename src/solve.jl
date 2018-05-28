function solvehook(m::JuMP.Model; suppress_warnings=false)
    data = getpolydata(m)
    for c in data.polyconstr
        c.delegate = getpolymodule(c).addpolyconstraint!(m, c.p, c.set, c.domain; c.kwargs...)
    end
    status = solve(m,suppress_warnings=suppress_warnings, ignore_solve_hook=true)
    return status
end
