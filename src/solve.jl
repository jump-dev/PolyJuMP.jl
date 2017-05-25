function solvehook(m::JuMP.Model; suppress_warnings=false)
    data = getpolydata(m)
    for c in data.polyconstr
        if c.nonnegative
            c.delegate = getpolymodule(c).addpolynonnegativeconstraint(m, c.p, c.domain)
        else
            c.delegate = getpolymodule(c).addpolyeqzeroconstraint(m, c.p, c.domain)
        end
    end
    status = solve(m,suppress_warnings=suppress_warnings, ignore_solve_hook=true)
    return status
end
