import StarAlgebras as SA
using Test, JuMP

function _model(optimizer::MOI.AbstractOptimizer)
    MOI.empty!(optimizer)
    return direct_model(optimizer)
end

function _model(factory)
    return Model(factory)
end

#"""
#    @test_suite setname subsets
#
#Defines a function `setname_test(model, config, exclude)` that runs the tests
#defined in the dictionary `setname_tests` with the model `model` and config
#`config` except the tests whose dictionary key is in `exclude`. If `subsets` is
#`true` then each test runs in fact multiple tests hence the `exclude` argument
#is passed as it can also contains test to be excluded from these subsets of
#tests.
#"""
macro test_suite(setname, subsets = false)
    testname = Symbol(string(setname) * "_test")
    testdict = Symbol(string(testname) * "s")
    if subsets
        runtest = :(f(model, config, exclude))
    else
        runtest = :(f(model, config))
    end
    return esc(
        :(
            function $testname(
                model, # could be ModelLike or an optimizer constructor
                config::$MOI.Test.Config,
                exclude::Vector{String} = String[],
            )
                for (name, f) in $testdict
                    if name in exclude
                        continue
                    end
                    @testset "$name" begin
                        $runtest
                    end
                end
            end
        ),
    )
end

function test_noc(model, F, S, n)
    @test MOI.get(model, MOI.NumberOfConstraints{F,S}()) == n
    @test length(MOI.get(model, MOI.ListOfConstraintIndices{F,S}())) == n
    @test ((F, S) in MOI.get(model, MOI.ListOfConstraintTypesPresent())) ==
          !iszero(n)
end

# Test deletion of bridge
function test_delete_bridge(
    model::Model,
    cref::ConstraintRef{Model,MOI.ConstraintIndex{F,S}},
    nvars::Int,
    nocs::Tuple,
) where {F,S}
    @test num_variables(model) == nvars
    @test length(all_variables(model)) == nvars
    test_noc(model, F, S, 1)
    for noc in nocs
        test_noc(model, noc...)
    end
    @test is_valid(model, cref)
    delete(model, cref)
    @test_throws MOI.InvalidIndex(index(cref)) delete(model, cref)
    @test !is_valid(model, cref)
    test_noc(model, F, S, 0)
    # As the bridge has been removed, if the constraints it has created where not removed, it wouldn't be there to decrease this counter anymore
    @test num_variables(model) == nvars
    @test length(all_variables(model)) == nvars
    for noc in nocs
        test_noc(model, noc...)
    end
end
