# PolyJuMP

| **PackageEvaluator** | **Build Status** |
|:--------------------:|:----------------:|
| [![][pkg-0.6-img]][pkg-0.6-url] | [![Build Status][build-img]][build-url] [![Build Status][winbuild-img]][winbuild-url] |
| [![][pkg-0.7-img]][pkg-0.7-url] | [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] |

PolyJuMP is an extension to [JuMP](https://github.com/JuliaOpt/JuMP.jl) for formulating and solving polynomial optimization problems. These problems can then be solved using [Sum of Squares Programming](https://github.com/JuliaOpt/SumOfSquares.jl).

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of [SumOfSquares](https://github.com/JuliaOpt/SumOfSquares.jl)' documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of [SumOfSquares](https://github.com/JuliaOpt/SumOfSquares.jl)' documentation.*

Some presentations on, or using, PolyJuMP:
  * Benoît Legat at the JuMP Meetup 2017 [[Slides](http://www.juliaopt.org/meetings/mit2017/legat.pdf)] [[Video](https://youtu.be/kyo72yWYr54)]
  * [Joey Huchette at SIAM Opt 2017](https://docs.google.com/presentation/d/1ASfjB1TdLJmYxT0b6rnyGh9eLbMc-66bTOt3_3yvc90/edit?usp=sharing)

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://juliaopt.github.io/SumOfSquares.jl/stable
[docs-latest-url]: https://juliaopt.github.io/SumOfSquares.jl/latest

[pkg-0.6-img]: http://pkg.julialang.org/badges/PolyJuMP_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=PolyJuMP
[pkg-0.7-img]: http://pkg.julialang.org/badges/PolyJuMP_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=PolyJuMP

[build-img]: https://travis-ci.org/JuliaOpt/PolyJuMP.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaOpt/PolyJuMP.jl
[winbuild-img]: https://ci.appveyor.com/api/projects/status/2y6dc0j2xk4aa4v7?svg=true
[winbuild-url]: https://ci.appveyor.com/project/JuliaOpt/polyjump-jl
[coveralls-img]: https://coveralls.io/repos/github/JuliaOpt/PolyJuMP.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaOpt/PolyJuMP.jl?branch=master
[codecov-img]: http://codecov.io/github/JuliaOpt/PolyJuMP.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaOpt/PolyJuMP.jl?branch=master
