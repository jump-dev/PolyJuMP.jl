_isequal(p, q) = all(JuMP.isequal_canonical.(coefficients(p), coefficients(q)))
function _isequal(x::AbstractArray, y::AbstractArray)
    size(x) == size(y) && all(_isequal.(x, y))
end

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

    function testcon(m, cref, S::Type, p, ineqs, eqs, basis=PolyJuMP.MonomialBasis, kwargs=[])
        @test cref isa JuMP.ConstraintRef{Model, <:MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}, <:S}}
        c = JuMP.constraint_object(cref)
        set = JuMP.moi_set(c)
        @test set isa S
        if set isa PolyJuMP.PlusMinusSet
            set = set.set
        end
        @test set.basis == basis
        if !isempty(kwargs)
            @test length(set.kwargs) == length(kwargs)
            for (i, kw) in enumerate(set.kwargs)
                @test kw.first == kwargs[i][1]
                @test kw.second == kwargs[i][2]
            end
        end
        if S == TestPolyModule.PosDefMatrix
            expected_p = set.y' * p * set.y
        else
            expected_p = p
        end
        # == between JuMP affine expression is not accurate, e.g. β + α != α + β
        # == 0 is not defined either
        # c.p and p can be matrices
        @test _isequal(polynomial(JuMP.jump_function(c), set.monomials),
                       expected_p)
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
    @testset "NonNeg" begin
        S = TestPolyModule.NonNeg
        testcon(m, @constraint(m, p >= q + 1, domain = @set y >= 1 && dom),
                S, p - q - 1, [y-1, x^3 + x*y^2 + y - 1], [x^2 + y^2 - 1])
        testcon(m, @constraint(m, p <= q),
                S, q - p, [], [])
        testcon(m, @constraint(m, q - p in PolyJuMP.NonNegPoly()),
                S, q - p, [], [])
        testcon(m, @constraint(m, p + q >= 0, domain = @set x == y^3),
                S, p + q, [], [x - y^3])
        @testset "Custom keyword" begin
            testcon(m, @constraint(m, p <= q, maxdegree=1),
                    S, q - p, [], [], MonomialBasis, [(:maxdegree, 1)])
        end
    end
    @testset "ZeroPolynomialSet" begin
        @testset "ZeroPolynomialSet{FullSpace}" begin
            S = PolyJuMP.ZeroPolynomialSet{FullSpace}
            testcon(m, @constraint(m, p == q),
                    S, p - q, [], [])
            testcon(m, @constraint(m, p - q in PolyJuMP.ZeroPoly()),
                    S, p - q, [], [])
        end
        S = PolyJuMP.ZeroPolynomialSet
        testcon(m, @constraint(m, p == q, domain = @set x == 1 && f(x, y)),
                S, p - q, [], [x - 1, x + y - 2])
        testcon(m, @constraint(m, p - q in PolyJuMP.ZeroPoly(), domain = @set x == 1 && f(x, y)),
                S, p - q, [], [x - 1, x + y - 2])
        S = PolyJuMP.PlusMinusSet
        testcon(m, @constraint(m, p == q, domain = dom),
                S, p - q, [x^3 + x*y^2 + y - 1],
                [x^2 + y^2 - 1])
    end
    @testset "PosDefMatrix" begin
        testcon(m, @SDconstraint(m, [p q; q 0] ⪰ [0 0; 0 p]),
                TestPolyModule.PosDefMatrix, [p q; q -p], [], [])
    end
end
