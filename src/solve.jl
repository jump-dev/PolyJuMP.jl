const MOIU = MOI.Utilities

function solvehook(m::JuMP.Model; suppress_warnings=false)
    data = getpolydata(m)
    for c in data.polyconstr
        c.delegate = getpolymodule(c).addpolyconstraint!(m, c.p, c.set, c.domain)
    end
    MOIU.resetsolver!(m, data.solver()) # temporary fix
    # TODO readd suppress_warnings when it is readded to JuMP
    status = solve(m, ignore_solve_hook=true)
    return status
end
