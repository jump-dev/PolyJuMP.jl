# PolyJuMP.jl

[![Build Status](https://github.com/jump-dev/PolyJuMP.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/PolyJuMP.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/PolyJuMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/PolyJuMP.jl)

[PolyJuMP.jl](https://github.com/jump-dev/PolyJuMP.jl) is a [JuMP](https://github.com/jump-dev/JuMP.jl)
extension for formulating and solving polynomial optimization problems.

These problems can then be solved using [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl).

## License

`PolyJuMP.jl` is licensed under the [MIT license](https://github.com/jump-dev/PolyJuMP.jl/blob/master/LICENSE.md).

## Installation

Install `PolyJuMP` using `Pkg.add`:
```julia
import Pkg
Pkg.add("PolyJuMP")
```

## Use with JuMP

To use QCQP solver with JuMP, use a nonconvex QCQP solver, e.g., `Gurobi.Optimizer` and `PolyJuMP.QCQP.Optimizer`:

```julia
using JuMP, PolyJuMP, Gurobi
model = Model(() -> PolyJuMP.QCQP.Optimizer(Gurobi.Optimizer))
```

To use KKT solver with JuMP, use solver of algebraic systems of equations implementing the [SemialgebraicSets interface](https://github.com/JuliaAlgebra/SemialgebraicSets.jl), e.g., `HomotopyContinuation.SemialgebraicSetsHCSolver` and `PolyJuMP.KKT.Optimizer`:

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
