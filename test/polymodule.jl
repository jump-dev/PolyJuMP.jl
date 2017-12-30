@testset "PolyModule" begin
    m = Model()
    # Triggers the creation of polydata
    @test isnull(PolyJuMP.getpolydata(m).polymodule)
    @test PolyJuMP.getpolymodule(m) == PolyJuMP.DefaultModule
    setpolymodule!(m, TestPolyModule)
    @test PolyJuMP.getpolymodule(m) == TestPolyModule
end
