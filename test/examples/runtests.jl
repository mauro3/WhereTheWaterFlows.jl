using Test

function _run_example_in_module(path::AbstractString)
    modname = Symbol("Example_", hash(path))
    m = Module(modname)
    Core.eval(m, :(using WhereTheWaterFlows))
    Base.include(m, path)
    return nothing
end

const EXAMPLES_DIR = normpath(joinpath(@__DIR__, "../../examples"))
const SMOKE_EXAMPLES = [
    joinpath(EXAMPLES_DIR, "wwf-simple.jl"),
    joinpath(EXAMPLES_DIR, "wwfs-simple.jl"),
    joinpath(EXAMPLES_DIR, "wwfr-simple.jl"),
]

@testset "examples smoke" begin
    @test !isempty(SMOKE_EXAMPLES)
    for fl in SMOKE_EXAMPLES
        @testset "$(relpath(fl, EXAMPLES_DIR))" begin
            @test_nowarn _run_example_in_module(fl)
        end
    end
end
