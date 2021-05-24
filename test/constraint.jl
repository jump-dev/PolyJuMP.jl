_isequal(p, q) = all(JuMP.isequal_canonical.(coefficients(p), coefficients(q)))
function _isequal(x::AbstractArray, y::AbstractArray)
    size(x) == size(y) && all(_isequal.(x, y))
end

# `MOI.Utilities.Model` canonicalizes the constraints so we need to
# canonicalize them as well for the printing tests.
function _canon(model, p::MP.APL)
    return MP.polynomial(map(MP.terms(p)) do t
        coef = MP.coefficient(t)
        moi = JuMP.moi_function(coef)
        jump = JuMP.jump_function(model, MOI.Utilities.canonical(moi))
        return MP.term(jump, MP.monomial(t))
    end)
end
_canon(model, p::Matrix) = _canon.(model, p)

@testset "@constraint macro with polynomials" begin
    m = Model()
    setpolymodule!(m, TestPolyModule)
    @variable m α
    @variable m β
    @polyvar x y
    p = α * x*y + β * x^2
    #q = MatPolynomial([α β; β α], [x])
    q = α*x^2 + β*x*y + α*y^2
    @test_macro_throws ErrorException @constraint(m, p)
    @test_macro_throws ErrorException @constraint(m, begin p >= 0 end)
    @test_macro_throws ErrorException @constraint(m, +(p, p, p))
    @test_macro_throws ErrorException @constraint(m, p >= 0, 1)
    #@test_macro_throws ErrorException @constraint(m, p >= 0, domain = (@set x >= -1 && x <= 1, domain = y >= -1 && y <= 1))
    @test_macro_throws ErrorException @constraint(m, p + 0, domain = (@set x >= -1 && x <= 1))

    function testcon(m, cref, S::Type, jump_set::PolyJuMP.PolynomialSet,
                     p, ineqs, eqs, basis=MB.MonomialBasis, kwargs=[])
        @test cref isa JuMP.ConstraintRef{
            Model, <:MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},
                                         <:S}}
        c = JuMP.constraint_object(cref)
        set = JuMP.moi_set(c)
        @test set isa S
        if set isa PolyJuMP.PlusMinusSet
            set = set.set
        end
        p_canon = _canon(m, p)
        expected_str = string(JuMP.function_string(REPLMode, p_canon), ' ',
                              JuMP._math_symbol(REPLMode, :in), ' ', jump_set)
        con = JuMP.constraint_object(cref)
        @test sprint(show, MIME"text/plain"(), cref) == expected_str
        expected_str = string("\$\$ ", JuMP.function_string(IJuliaMode, p_canon), ' ',
                              JuMP._math_symbol(IJuliaMode, :in), ' ', jump_set,
                              " \$\$")
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
        @test _isequal(JuMP.reshape_vector(JuMP.jump_function(c),
                                           JuMP.shape(c)), p)
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

    f(x, y) = @set x + y == 2
    dom = @set x^2 + y^2 == 1 && x^3 + x*y^2 + y >= 1

    @testset "Printing" begin
        in_sym = JuMP._math_symbol(REPLMode, :in)
        eqref = @constraint(m, p == q)
        @test sprint(show, MIME"text/plain"(), eqref) == "(-α + β)x² + (α - β)xy + (-α)y² $in_sym PolyJuMP.ZeroPoly()"
        @test sprint(show, MIME"text/latex"(), eqref) == "\$\$ (-α + β)x^{2} + (α - β)xy + (-α)y^{2} \\in PolyJuMP.ZeroPoly() \$\$"
        sdref = @constraint(m, [p q; q p] in PSDCone())
        if VERSION < v"1.4-"
            @test sprint(show, MIME"text/plain"(), sdref) == "[(β)x² + (α)xy          (α)x² + (β)xy + (α)y²;\n (α)x² + (β)xy + (α)y²  (β)x² + (α)xy        ] $in_sym Main.TestPolyModule.TestPosDefMatrix()"
        else
            @test sprint(show, MIME"text/plain"(), sdref) == "[(β)x² + (α)xy          (α)x² + (β)xy + (α)y²;\n (α)x² + (β)xy + (α)y²  (β)x² + (α)xy] $in_sym Main.TestPolyModule.TestPosDefMatrix()"
        end
        @test sprint(show, MIME"text/latex"(), sdref) == "\$\$ \\begin{bmatrix}\n(β)x^{2} + (α)xy & (α)x^{2} + (β)xy + (α)y^{2}\\\\\n(α)x^{2} + (β)xy + (α)y^{2} & (β)x^{2} + (α)xy\\\\\n\\end{bmatrix} \\in Main.TestPolyModule.TestPosDefMatrix() \$\$"
    end

    @testset "NonNeg" begin
        S = TestPolyModule.NonNeg
        jump_set = TestPolyModule.TestNonNeg()
        testcon(m, @constraint(m, p >= q + 1, domain = @set y >= 1 && dom),
                S, jump_set, p - q - 1, [y - 1, x^3 + x*y^2 + y - 1],
                [x^2 + y^2 - 1])
        testcon(m, @constraint(m, p <= q),
                S, jump_set, -p + q, [], [])
        testcon(m, @constraint(m, q - p in PolyJuMP.NonNegPoly()),
                S, jump_set, q - p, [], [])
        testcon(m, @constraint(m, p + q >= 0, domain = @set x == y^3),
                S, jump_set, p + q, [], [x - y^3])
        @testset "Custom keyword" begin
            testcon(m, @constraint(m, p <= q, maxdegree=1),
                    S, jump_set, -p + q, [], [], MB.MonomialBasis,
                    [(:maxdegree, 1)])
        end
    end
    @testset "ZeroPolynomialSet" begin
        jump_set = PolyJuMP.ZeroPoly()
        @testset "ZeroPolynomialSet{FullSpace}" begin
            S = PolyJuMP.ZeroPolynomialSet{FullSpace}
            testcon(m, @constraint(m, p == q),
                    S, jump_set, p - q, [], [])
            testcon(m, @constraint(m, p - q in PolyJuMP.ZeroPoly()),
                    S, jump_set, p - q, [], [])
        end
        S = PolyJuMP.ZeroPolynomialSet
        testcon(m, @constraint(m, p == q, domain = @set x == 1 && f(x, y)),
                S, jump_set, p - q, [], [x + y - 2, x - 1])
        testcon(m, @constraint(m, p - q in PolyJuMP.ZeroPoly(), domain = @set x == 1 && f(x, y)),
                S, jump_set, p - q, [], [x + y - 2, x - 1])
        S = PolyJuMP.PlusMinusSet
        testcon(m, @constraint(m, p == q, domain = dom),
                S, jump_set, p - q, [x^3 + x*y^2 + y - 1],
                [x^2 + y^2 - 1])
    end
    @testset "PosDefMatrix" begin
        jump_set = TestPolyModule.TestPosDefMatrix()
        testcon(m, @SDconstraint(m, [p q; q 0] ⪰ [0 0; 0 p]),
                TestPolyModule.PosDefMatrix, jump_set, [p q; q -p], [], [])
    end
end
