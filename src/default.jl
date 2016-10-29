const polymodules = [:SumOfSquares]

function getdefaultpolymodule()
    for pkgname in polymodules
        try
            eval(Expr(:import,pkgname))
        catch
            if isdir(Pkg.dir((string(pkgname))))
                warn("Package ",string(pkgname),
                " is installed but couldn't be loaded. ",
                "You may need to run `Pkg.build(\"$pkgname\")`")
            end
            continue
        end
        return eval(pkgname)
    end
    return nothing
end
