using SumOfSquares

const polyhedra_test = joinpath(Pkg.dir("SumOfSquares"), "test")

sostest(s) = include(joinpath(polyhedra_test, "$s.jl"))

sostest("solvers")

sostest("certificate")

sostest("motzkin")

# SOSTools demos
sostest("sospoly")
sostest("sosdemo2")
sostest("sosdemo3")
sostest("sosdemo4")
sostest("sosdemo5")
sostest("sosdemo6")
sostest("domain")
sostest("sosmatrix")
