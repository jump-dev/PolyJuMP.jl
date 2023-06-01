module TestConstraint

using Test

using MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases
using SemialgebraicSets

using JuMP
using PolyJuMP

include("utilities.jl")
include("testpolymodule.jl")

_isequal(p, q) = all(JuMP.isequal_canonical.(coefficients(p), coefficients(q)))
function _isequal(x::AbstractArray, y::AbstractArray)
    return size(x) == size(y) && all(_isequal.(x, y))
end

_con_constant(a::Real) = convert(GenericAffExpr{Float64,VariableRef}, a)
_con_constant(a::Complex) = convert(GenericAffExpr{ComplexF64,VariableRef}, a)
_con_constant(a) = a

# `MOI.Utilities.Model` canonicalizes the constraints so we need to
# canonicalize them as well for the printing tests.
function _canon(model, p::MP.AbstractPolynomialLike)
    return MP.polynomial(
        map(MP.terms(p)) do t
            coef = _con_constant(MP.coefficient(t))
            moi = JuMP.moi_function(coef)
            jump = JuMP.jump_function(model, MOI.Utilities.canonical(moi))
            return MP.term(jump, MP.monomial(t))
        end,
    )
end
_canon(model, p::Matrix) = _canon.(model, p)

function _test_constraint(
    m,
    cref,
    S::Type,
    jump_set::PolyJuMP.PolynomialSet,
    p,
    ineqs,
    eqs,
    basis = MB.MonomialBasis,
    kwargs = [];
    T = Float64,
)
    @test cref isa JuMP.ConstraintRef{
        Model,
        <:MOI.ConstraintIndex{MOI.VectorAffineFunction{T},<:S},
    }
    c = JuMP.constraint_object(cref)
    set = JuMP.moi_set(c)
    @test set isa S
    if set isa PolyJuMP.PlusMinusSet
        set = set.set
    end
    p_canon = _canon(m, p)
    in_str = Sys.iswindows() ? "in" : "∈"
    expected_str = string(
        JuMP.function_string(MIME"text/plain"(), p_canon),
        ' ',
        in_str,
        ' ',
        jump_set,
    )
    con = JuMP.constraint_object(cref)
    @test sprint(show, MIME"text/plain"(), cref) == expected_str
    expected_str = string(
        "\$\$ ",
        JuMP.function_string(MIME"text/latex"(), p_canon),
        ' ',
        "\\in ",
        jump_set,
        " \$\$",
    )
    @test sprint(show, MIME"text/latex"(), cref) == expected_str
    @test set.basis == basis
    if !isempty(kwargs)
        @test length(set.kwargs) == length(kwargs)
        for (i, kw) in enumerate(set.kwargs)
            @test kw.first == kwargs[i][1]
            @test kw.second == kwargs[i][2]
        end
    end
    # == between JuMP affine expression is not accurate, e.g. β + α != α + β
    # == 0 is not defined either
    # c.p and p can be matrices
    @test _isequal(
        JuMP.reshape_vector(JuMP.jump_function(c), JuMP.shape(c)),
        p_canon,
    )
    if isempty(ineqs)
        if isempty(eqs)
            @test isa(set.domain, FullSpace)
        else
            @test isa(set.domain, AlgebraicSet)
            @test equalities(set.domain) == eqs
        end
    else
        @test isa(set.domain, BasicSemialgebraicSet)
        @test inequalities(set.domain) == ineqs
        @test equalities(set.domain) == eqs
    end
end

function test_errors(var)
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
    m = Model()
    setpolymodule!(m, DummyPolyModule)
    @variable m α
    @variable m β
    p = α * x * y + β * x^2
    q = α * x^2 + β * x * y + α * y^2
    @test_macro_throws ErrorException @constraint(m, p)
    @test_macro_throws ErrorException @constraint(m, begin
        p >= 0
    end)
    @test_macro_throws ErrorException @constraint(m, +(p, p, p))
    @test_macro_throws ErrorException @constraint(m, p >= 0, 1)
    #@test_macro_throws ErrorException @constraint(m, p >= 0, domain = (@set x >= -1 && x <= 1, domain = y >= -1 && y <= 1))
    @test_macro_throws ErrorException @constraint(
        m,
        p + 0,
        domain = (@set x >= -1 && x <= 1)
    )
end

function test_printing(var)
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
    m = Model()
    setpolymodule!(m, DummyPolyModule)
    @variable m α
    @variable m β
    p = α * x * y + β * x^2
    q = α * x^2 + β * x * y + α * y^2

    in_sym = Sys.iswindows() ? "in" : "∈"
    eqref = @constraint(m, p == q)
    @test sprint(show, MIME"text/plain"(), eqref) ==
          "(-α)y² + (α - β)xy + (-α + β)x² $in_sym PolyJuMP.ZeroPoly()"
    @test sprint(show, MIME"text/latex"(), eqref) ==
          "\$\$ (-α)y^{2} + (α - β)xy + (-α + β)x^{2} \\in PolyJuMP.ZeroPoly() \$\$"
    sdref = @constraint(m, [p q; q p] in PSDCone())
    @test sprint(show, MIME"text/plain"(), sdref) ==
          "[(α)xy + (β)x²          (α)y² + (β)xy + (α)x²;\n (α)y² + (β)xy + (α)x²  (α)xy + (β)x²] $in_sym $DummyPolyModule.DummyPosDefMatrix()"
    @test sprint(show, MIME"text/latex"(), sdref) ==
          "\$\$ \\begin{bmatrix}\n(α)xy + (β)x^{2} & (α)y^{2} + (β)xy + (α)x^{2}\\\\\n(α)y^{2} + (β)xy + (α)x^{2} & (α)xy + (β)x^{2}\\\\\n\\end{bmatrix} \\in $DummyPolyModule.DummyPosDefMatrix() \$\$"
end

function test_NonNeg(var)
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
    m = Model()
    setpolymodule!(m, DummyPolyModule)
    @variable m α
    @variable m β
    p = α * x * y + β * x^2
    q = α * x^2 + β * x * y + α * y^2
    f(x, y) = @set x + y == 2
    dom = @set x^2 + y^2 == 1 && x^3 + x * y^2 + y >= 1
    S = DummyPolyModule.NonNeg
    jump_set = DummyPolyModule.DummyNonNeg()
    _test_constraint(m, @constraint(m, x >= y), S, jump_set, x - y, [], [])
    _test_constraint(
        m,
        @constraint(m, im * x >= y),
        S,
        jump_set,
        im * x - y,
        [],
        [],
        T = ComplexF64,
    )
    _test_constraint(
        m,
        @constraint(m, p >= q + 1, domain = @set y >= 1 && dom),
        S,
        jump_set,
        p - q - 1,
        [y - 1, x^3 + x * y^2 + y - 1],
        [x^2 + y^2 - 1],
    )
    _test_constraint(m, @constraint(m, p <= q), S, jump_set, -p + q, [], [])
    _test_constraint(
        m,
        @constraint(m, q - p in PolyJuMP.NonNegPoly()),
        S,
        jump_set,
        q - p,
        [],
        [],
    )
    _test_constraint(
        m,
        @constraint(m, p + q >= 0, domain = @set x == y^3),
        S,
        jump_set,
        p + q,
        [],
        [x - y^3],
    )
    @testset "Custom keyword" begin
        _test_constraint(
            m,
            @constraint(m, p <= q, maxdegree = 1),
            S,
            jump_set,
            -p + q,
            [],
            [],
            MB.MonomialBasis,
            [(:maxdegree, 1)],
        )
    end
end

function test_ZeroPolynomialSet(var)
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
    m = Model()
    setpolymodule!(m, DummyPolyModule)
    @variable m α
    @variable m β
    p = α * x * y + β * x^2
    q = α * x^2 + β * x * y + α * y^2
    f(x, y) = @set x + y == 2
    dom = @set x^2 + y^2 == 1 && x^3 + x * y^2 + y >= 1
    jump_set = PolyJuMP.ZeroPoly()
    @testset "ZeroPolynomialSet{FullSpace}" begin
        S = PolyJuMP.ZeroPolynomialSet{FullSpace}
        _test_constraint(m, @constraint(m, p == q), S, jump_set, p - q, [], [])
        @test PolyJuMP.Bridges.Constraint.ZeroPolynomialBridge in m.bridge_types
        _test_constraint(
            m,
            @constraint(m, p - q in PolyJuMP.ZeroPoly()),
            S,
            jump_set,
            p - q,
            [],
            [],
        )
        _test_constraint(m, @constraint(m, x == y), S, jump_set, x - y, [], [])
        _test_constraint(
            m,
            @constraint(m, x == im * y),
            S,
            jump_set,
            x - im * y,
            [],
            [],
            T = ComplexF64,
        )
    end
    S = PolyJuMP.ZeroPolynomialSet
    _test_constraint(
        m,
        @constraint(m, p == q, domain = @set x == 1 && f(x, y)),
        S,
        jump_set,
        p - q,
        [],
        [x + y - 2, x - 1],
    )
    @test PolyJuMP.Bridges.Constraint.ZeroPolynomialInAlgebraicSetBridge in
          m.bridge_types
    _test_constraint(
        m,
        @constraint(
            m,
            p - q in PolyJuMP.ZeroPoly(),
            domain = @set x == 1 && f(x, y)
        ),
        S,
        jump_set,
        p - q,
        [],
        [x + y - 2, x - 1],
    )
    S = PolyJuMP.PlusMinusSet
    return _test_constraint(
        m,
        @constraint(m, p == q, domain = dom),
        S,
        jump_set,
        p - q,
        [x^3 + x * y^2 + y - 1],
        [x^2 + y^2 - 1],
    )
end

function test_PosDefMatrix(var)
    m = Model()
    setpolymodule!(m, DummyPolyModule)
    @variable m α
    @variable m β
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
    p = α * x * y + β * x^2
    q = α * x^2 + β * x * y + α * y^2

    jump_set = DummyPolyModule.DummyPosDefMatrix()
    _test_constraint(
        m,
        @constraint(m, [p q; q 0] - [0 0; 0 p] in PSDCone()),
        DummyPolyModule.PosDefMatrix,
        jump_set,
        [p q; q -p],
        [],
        [],
    )
    _test_constraint(
        m,
        @constraint(m, [p q; q 0] >= [0 0; 0 p], PSDCone()),
        DummyPolyModule.PosDefMatrix,
        jump_set,
        [p q; q -p],
        [],
        [],
    )
    return _test_constraint(
        m,
        @constraint(m, [0 0; 0 p] <= [p q; q 0], PSDCone()),
        DummyPolyModule.PosDefMatrix,
        jump_set,
        [p q; q -p],
        [],
        [],
    )
end

function runtests(var)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)(var)
            end
        end
    end
end

end

import DynamicPolynomials
DynamicPolynomials.@polyvar(x)
TestConstraint.runtests(x)
import TypedPolynomials
TypedPolynomials.@polyvar(y)
TestConstraint.runtests(y)
