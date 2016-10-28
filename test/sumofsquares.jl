using SumOfSquares

const polyhedra_test = joinpath(Pkg.dir("SumOfSquares"), "test")

sostest(s) = include(joinpath(polyhedra_test, "$s.jl"))

sostest("solvers")

sostest("motzkin")
sostest("choi")

# SOSTools demos
sostest("sosdemo1")
sostest("sosdemo2")
sostest("sosdemo3")
sostest("sosdemo4")
#sostest("sosdemo5")
#sostest("sosdemo6")
sostest("sosdemo9")
sostest("sosdemo10")
