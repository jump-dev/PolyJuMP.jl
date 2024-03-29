# PolyJuMP.jl

[![Build Status](https://github.com/jump-dev/PolyJuMP.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/PolyJuMP.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/PolyJuMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/PolyJuMP.jl)

[PolyJuMP.jl](https://github.com/jump-dev/PolyJuMP.jl) is a [JuMP](https://github.com/jump-dev/JuMP.jl)
extension for formulating and solving polynomial optimization problems.
This extension includes the following:

* Polynomial functions on JuMP decisions variables. These can be solved with the `PolyJuMP.QCQP.Optimizer` or `PolyJuMP.KKT.Optimizer`.
* Constraints that a polynomial is nonnegative where the coefficients of the polynomials depend on JuMP decision variables.
  These nonnegativity constraints can be reformulated using sufficient conditions using `PolyJuMP.RelativeEntropy` submodule or [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl).

## License

`PolyJuMP.jl` is licensed under the [MIT license](https://github.com/jump-dev/PolyJuMP.jl/blob/master/LICENSE.md).

## Installation

Install `PolyJuMP` using `Pkg.add`:
```julia
import Pkg
Pkg.add("PolyJuMP")
```

## Use with JuMP

### Polynomial nonnegativity constraints

`PolyJuMP` allows encoding a constraint that a polynomial should be nonnegative for all values of some
symbolic variables defined with `DynamicPolynomials.@polyvar` or `TypedPolynomials.@polyvar` as follows.
For instance, the following constrains the JuMP decision variable `a` to be such that
`a * x * y^2 + y^3 - a * x` is nonnegative for all real values of `x` and `y`:
```julia
using DynamicPolynomials
@polyvar x y
using JuMP
model = Model()
@variable(model, a)
@constraint(model, a * x * y^2 + y^3 >= a * x)
```
Determining the nonnegativity of a multivariate polynomial is however NP-hard so sufficient conditions
are used instead.
You need to specify which sufficient condition is used explicitly.
To use Sum-of-Arithmetic-Geometric-Exponentials (SAGE), use
```julia
import PolyJuMP
PolyJuMP.setpolymodule!(model, PolyJuMP.SAGE)
```
To use Sum-of-Squares (SOS), use
```julia
import SumOfSquares
PolyJuMP.setpolymodule!(model, SumOfSquares)
```
or replace `model = Model()` by `model = SOSModel()`.

Alternatively, the nonnegativity constraint can be explicit:
```julia
@constraint(model, a * x * y^2 + y^3 - a * x in PolyJuMP.SAGE.Polynomials())
@constraint(model, a * x * y^2 + y^3 - a * x in SumOfSquares.SOSCone())
```
This allows mixing SAGE and SOS constraints in the same model.

### Polynomial optimization

PolyJuMP also allows solving polynomial optimization problems using the `QCQP` and `KKT` solvers.
Polynomial optimization problems do not involve any symbolic variables from DynamicPolynomials or TypedPolynomials,
instead all variables are JuMP decision variables.

The `QCQP` solver is parametrized by a nonconvex QCQP inner solver.
It reformulates the polynomial optimization problem into a nonconvex `QCQP`
and relies on the inner solver to solve it.
For instance, to use the `QCQP` solver with JuMP with `Gurobi.Optimizer` as inner solver, use:
```julia
using JuMP, PolyJuMP, Gurobi
model = Model(() -> PolyJuMP.QCQP.Optimizer(Gurobi.Optimizer))
```

The `KKT` solver is parametrized by an inner solver of algebraic systems of equations implementing the [SemialgebraicSets interface](https://github.com/JuliaAlgebra/SemialgebraicSets.jl).
It reformulates the polynomial optimization problem into a system of polynomial equations
and relies on the inner solver to solve it.
For instance, to use the `QCQP` solver with JuMP with
`HomotopyContinuation.SemialgebraicSetsHCSolver` as inner solver, use:
```julia
using JuMP, PolyJuMP, HomotopyContinuation
model = Model(optimizer_with_attributes(
    PolyJuMP.KKT.Optimizer,
    "solver" => HomotopyContinuation.SemialgebraicSetsHCSolver(),
))
```

## Documentation

Documentation for `PolyJuMP.jl` is included in the
[documentation for SumOfSquares.jl](https://jump.dev/SumOfSquares.jl/stable).
