setpolymodule!(m::JuMP.Model, pm::Module) = setpolymodule!(getpolydata(m), pm)
setpolymodule!(data::PolyData, pm::Module) = data.polymodule = pm

getpolymodule(m::JuMP.Model) = getpolymodule(getpolydata(m))
function getpolymodule(data::PolyData)
    if isnull(data.polymodule)
        return PolyJuMP
    end
    get(data.polymodule)
end
