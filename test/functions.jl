@testset "Functions" begin
    @polyvar x
    α = MOI.VariableIndex(1)
    f = MOI.SingleVariable(α)
    @test f * x isa AbstractTerm{MOI.SingleVariable}
    @test x * f isa AbstractTerm{MOI.SingleVariable}
    @test f * x^2 isa AbstractTerm{MOI.SingleVariable}
    @test x^2 * f isa AbstractTerm{MOI.SingleVariable}
    @test f * (1x) isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test (1x) * f isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test f * (1.0x) isa AbstractTerm{MOI.ScalarAffineFunction{Float64}}
    @test (1.0x) * f isa AbstractTerm{MOI.ScalarAffineFunction{Float64}}
    @test f * (1x + 2) isa AbstractPolynomial{MOI.ScalarAffineFunction{Int}}
    @test (1x + 2) * f isa AbstractPolynomial{MOI.ScalarAffineFunction{Int}}
    @test f * (1.0x + 2.0) isa AbstractPolynomial{MOI.ScalarAffineFunction{Float64}}
    @test (1.0x + 2.0) * f isa AbstractPolynomial{MOI.ScalarAffineFunction{Float64}}
end
