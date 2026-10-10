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

PolyJuMP also allows solving polynomial optimization problems using the `QCQP`, `KKT` and `SAGE` solvers.
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

The `SAGE` solver computes a bound on the optimal objective value using the
SAGE relaxation of the problem [CP16, MCW21]: it certifies the nonnegativity
of the Lagrangian with the SAGE cone, using one SAGE multiplier for each
inequality constraint and one free polynomial multiplier for each equality
constraint.
It is parametrized by an inner solver for the resulting relative entropy
program. For instance, to compute a lower bound on the minimum of the Motzkin
polynomial with `ECOS.Optimizer` as inner solver, use:
```julia
using JuMP, PolyJuMP, ECOS
model = Model(() -> PolyJuMP.SAGE.Optimizer(ECOS.Optimizer))
@variable(model, x)
@variable(model, y)
@objective(model, Min, x^4 * y^2 + x^2 * y^4 + 1 - 3 * x^2 * y^2)
optimize!(model)
objective_bound(model) # ≈ 0
```
The bound is returned as `MOI.ObjectiveBound`. In addition, candidate
solutions are recovered from the dual of the SAGE constraint, which is a
vector of pseudo-moments, following [MCW21, Section 4.2] (see also its
reference implementation `poly_solrec` in
[sageopt](https://github.com/rileyjmurray/sageopt)); `result_count(model)`
gives the number of candidates found, sorted by feasibility and objective
value. In the example above, the four minimizers `(±1, ±1)` are recovered:
```julia
value(x; result = 1), value(y; result = 1) # ≈ (1, 1)
```
The maximum degree of the multiplier of a constraint is chosen with the
`PolyJuMP.MultiplierMaxdegree` constraint attribute:
```julia
@constraint(model, con, x^2 >= 1)
MOI.set(model, PolyJuMP.MultiplierMaxdegree(), con, 2)
```
The `SumOfSquares.Optimizer` of [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl)
is the analogous solver certifying the nonnegativity of the Lagrangian with
the SOS cone instead; increasing the `PolyJuMP.MultiplierMaxdegree` attributes
then gives the higher levels of the Lasserre hierarchy and the candidate
solutions are recovered from the atoms of the moment matrix given by the dual
of the SOS constraint, using
[MultivariateMoments](https://github.com/JuliaAlgebra/MultivariateMoments.jl).

[CP16] Chandrasekaran, Venkat, and Parikshit Shah.
*Relative entropy relaxations for signomial optimization.*
SIAM Journal on Optimization 26.2 (2016): 1147-1173.

[MCW21] Murray, Riley, Venkat Chandrasekaran, and Adam Wierman.
*Signomials and polynomial optimization via relative entropy and partial
dualization.* Mathematical Programming Computation 13 (2021): 257-295.
[arXiv:1907.00814](https://arxiv.org/abs/1907.00814)

## Documentation

Documentation for `PolyJuMP.jl` is included in the
[documentation for SumOfSquares.jl](https://jump.dev/SumOfSquares.jl/stable).
