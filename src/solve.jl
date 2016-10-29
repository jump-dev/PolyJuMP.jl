function solvehook(m::JuMP.Model; suppress_warnings=false)
    data = getpolydata(m)
    for c in data.polyconstr
        if c.nonnegative
            c.delegate = c.polymodule.addpolynonnegativeconstraint(m, c.p)
        else
            c.delegate = c.polymodule.addpolyeqzeroconstraint(m, c.p)
        end
    end
    status = solve(m,suppress_warnings=suppress_warnings, ignore_solve_hook=true)
    return status
end
