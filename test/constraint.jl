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

    function testcon(m, cref, set::ZeroPoly, dom)
        @test string(cref) == "PolyJuMP constraint"
        @test isa(cref, JuMP.ConstraintRef{Model, <:Union{PolyJuMP.ZeroConstraint, PolyJuMP.ZeroConstraintWithDomain}})
        delegate = PolyJuMP.getdelegate(cref)
        # TODO test VectorAffineFunction vs VectorQuadraticFunction
        #@test delegate.x == x # TODO
        if dom
            @test delegate isa PolyJuMP.ZeroConstraintWithDomain
        else
            @test delegate isa PolyJuMP.ZeroConstraint
        end
    end
    function testcon(m, cref, set, p, ineqs, eqs, basis=PolyJuMP.MonomialBasis, kwargs=[])
        @test string(cref) == "PolyJuMP constraint"
        @test isa(cref, JuMP.ConstraintRef{Model, TestPolyModule.TestConstraint})
        c = PolyJuMP.getdelegate(cref)
        @test c.basis == basis
        @test c.set == set
        @test length(c.kwargs) == length(kwargs)
        for (i, kw) in enumerate(c.kwargs)
            @test kw.first == kwargs[i][1]
            @test kw.second == kwargs[i][2]
        end
        # == between JuMP affine expression is not accurate, e.g. β + α != α + β
        # == 0 is not defined either
        # c.p and p can be matrices
        @test _isequal(c.p, p)
        if isempty(ineqs)
            if isempty(eqs)
                @test isa(c.domain, FullSpace)
            else
                @test isa(c.domain, AlgebraicSet)
                @test equalities(c.domain) == eqs
            end
        else
            @test isa(c.domain, BasicSemialgebraicSet)
            @test inequalities(c.domain) == ineqs
            @test equalities(c.domain) == eqs
        end
    end

    f(x, y) = @set x + y == 2
    dom = @set x^2 + y^2 == 1 && x^3 + x*y^2 + y >= 1
    testcon(m, @constraint(m, p >= q + 1, domain = @set y >= 1 && dom), TestPolyModule.TestNonNegConstraint(), p - q - 1, [y-1, x^3 + x*y^2 + y - 1], [x^2 + y^2 - 1])
    testcon(m, @constraint(m, p <= q), TestPolyModule.TestNonNegConstraint(), q - p, [], [])
    testcon(m, @constraint(m, q - p in NonNegPoly()), TestPolyModule.TestNonNegConstraint(), q - p, [], [])
    testcon(m, @constraint(m, p + q >= 0, domain = @set x == y^3), TestPolyModule.TestNonNegConstraint(), p + q, [], [x - y^3])
    testcon(m, @constraint(m, p == q, domain = @set x == 1 && f(x, y)), ZeroPoly(), false)
    testcon(m, @constraint(m, p == q, domain = dom), ZeroPoly(), true)
    testcon(m, @constraint(m, p - q in ZeroPoly(), domain = @set x == 1 && f(x, y)), ZeroPoly(), false)
    testcon(m, @SDconstraint(m, [p q; q 0] ⪰ [0 0; 0 p]), TestPolyModule.TestNonNegMatrixConstraint(), [p q; q -p], [], [])
    testcon(m, @constraint(m, p <= q, maxdegree=1), TestPolyModule.TestNonNegConstraint(), q - p, [], [], MonomialBasis, [(:maxdegree, 1)])
end
