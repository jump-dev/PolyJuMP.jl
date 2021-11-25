@testset "Functions" begin
    @polyvar x
    α = MOI.VariableIndex(1)
    @test α * x isa AbstractTerm{MOI.VariableIndex}
    @test x * α isa AbstractTerm{MOI.VariableIndex}
    @test α * x^2 isa AbstractTerm{MOI.VariableIndex}
    @test x^2 * α isa AbstractTerm{MOI.VariableIndex}
    @test (1α) * x^2 isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test (1x)^2 * α isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test α * (1x) isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test (1x) * α isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test α * (1.0x) isa AbstractTerm{MOI.ScalarAffineFunction{Float64}}
    @test (1.0x) * α isa AbstractTerm{MOI.ScalarAffineFunction{Float64}}
    @test α * (1x + 2) isa AbstractPolynomial{MOI.ScalarAffineFunction{Int}}
    @test (1x + 2) * α isa AbstractPolynomial{MOI.ScalarAffineFunction{Int}}
    @test α * (1.0x + 2.0) isa AbstractPolynomial{MOI.ScalarAffineFunction{Float64}}
    @test (1.0x + 2.0) * α isa AbstractPolynomial{MOI.ScalarAffineFunction{Float64}}
end
