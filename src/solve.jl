function solvehook(m::JuMP.Model; suppress_warnings=false)
    for c in getpolyconstr(m)
        c.delegate = addpolyconstraint!(m, c.p, c.set, c.domain; c.kwargs...)
    end
    solve(m, suppress_warnings=suppress_warnings, ignore_solve_hook=true)
end
