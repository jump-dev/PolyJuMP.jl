using MultivariatePolynomials
using JuMP
using PolyJuMP
using Base.Test

if isdir(Pkg.dir("SumOfSquares"))
    include("sumofsquares.jl")
end
