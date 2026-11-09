module GLOrderingTests
using Test, WhereTheWaterFlows
const WWF = WhereTheWaterFlows
const WWFS = WWF.Subglacially
const CI = CartesianIndex
include("gl-order-fixture.jl")
const original_fixture = gl_order_fixture()

is_d4(path) = all(sum(abs, Tuple(J - I)) == 1 for (I, J) in zip(path, path[2:end]))

# Independent synthetic boundary definition, including diagonal exterior contacts.
function boundary(mask)
    [I for I in CartesianIndices(mask) if mask[I] &&
        any(d -> checkbounds(Bool, mask, I + d) && !mask[I + d], WWFS.D8_POS)]
end

@testset "Original GL fixture, including every starting cell" begin
    gl, mask = original_fixture.gl, original_fixture.mask
    g = Set(gl)
    inside = WWFS.walk_gl_inside(gl, mask)
    @test length(gl) == 100
    @test length(inside) == 88
    @test allunique(inside)
    @test is_d4(inside)
    @test CI(22, 28) ∉ inside
    @test WWFS.step_gl(CI(21, 27), CI(21, 28), g, mask) == CI(21, 29)
    @test WWFS.step_gl(CI(21, 26), CI(21, 27), g, mask;
                      region=:outside, side=:right) == CI(22, 28)
    for start in gl
        outside = WWFS.walk_gl_outside(gl, mask, start)
        @test first(outside) == start
        @test allunique(outside)
        @test Set(outside) == g
        result = WWFS.walk_gl(gl, mask, start)
        @test Set(result.inside) == Set(inside)
        @test allunique(result.inside)
        @test is_d4(result.inside)
        @test result.outside == outside
        unculled = WWFS.walk_gl(gl, mask, start; cull_inside=false)
        @test unculled.inside == [result.inside]
        @test unculled.outside == outside
        if start in inside
            path = WWFS.walk_gl_inside(gl, mask, start)
            @test Set(path) == Set(inside)
            @test allunique(path)
            @test is_d4(path)
        else
            @test_throws ErrorException WWFS.walk_gl_inside(gl, mask, start)
        end
    end
    @test_throws ErrorException WWFS.walk_gl_outside(gl, mask, CI(1, 1))
end

@testset "GL reflection, transpose, and input order" begin
    for reflect in (false, true), transpose in (false, true)
        mask = reflect ? reverse(original_fixture.mask; dims=1) : original_fixture.mask
        gl = reflect ? [CI(size(mask, 1) + 1 - I[1], I[2]) for I in original_fixture.gl] :
                       original_fixture.gl
        if transpose
            mask = permutedims(mask)
            gl = [CI(I[2], I[1]) for I in gl]
        end
        for cells in (gl, reverse(gl))
            inside = WWFS.walk_gl_inside(cells, mask)
            outside = WWFS.walk_gl_outside(cells, mask)
            @test length(inside) == 88
            @test allunique(inside)
            @test is_d4(inside)
            @test allunique(outside)
            @test Set(outside) == Set(gl)
        end
    end
end

@testset "Closed GL loop" begin
    mask = falses(7, 7)
    mask[2:6, 2:6] .= true
    gl = boundary(mask)
    for start in gl
        inside = WWFS.walk_gl_inside(gl, mask, start)
        @test Set(inside) == Set(gl)
        @test allunique(inside)
        @test is_d4([inside; first(inside)])
        @test Set(WWFS.walk_gl_outside(gl, mask, start)) == Set(gl)
    end
end

@testset "Long narrow GL spur requires outside backtracking" begin
    mask = falses(12, 13)
    mask[2:6, 2:12] .= true
    mask[7:10, 7] .= true
    gl = boundary(mask)
    inside = WWFS.walk_gl_inside(gl, mask)
    @test is_d4(inside)
    @test allunique(inside)
    @test all(I[1] <= 6 for I in inside)
    for start in (first(gl), CI(10, 7))
        outside = WWFS.walk_gl_outside(gl, mask, start)
        @test allunique(outside)
        @test Set(outside) == Set(gl)
    end
    outside = WWFS.walk_gl_outside(gl, mask)
    @test any(maximum(abs, Tuple(J - I)) > 1
              for (I, J) in zip(outside, outside[2:end]))
end

@testset "Diagonal GL attachment with a separate interior" begin
    mask = falses(11, 11)
    mask[2:6, 2:6] .= true
    mask[7:9, 7:9] .= true
    gl = boundary(mask)
    for (start, n) in ((CI(2, 2), 16), (CI(7, 7), 8))
        inside = WWFS.walk_gl_inside(gl, mask, start)
        @test length(inside) == n
        @test allunique(inside)
        @test is_d4(inside)
        outside = WWFS.walk_gl_outside(gl, mask, start)
        @test allunique(outside)
        @test Set(outside) == Set(gl)
        result = WWFS.walk_gl(gl, mask, start)
        @test length(result.inside) == 16
        @test allunique(result.inside)
        @test is_d4([result.inside; first(result.inside)])
        @test result.outside == outside
    end
end

@testset "Added island: longest inside GL walk from every start" begin
    mask = copy(original_fixture.mask)
    mask[22:35, 40:46] .= true
    gl = boundary(mask)
    expected = WWFS.walk_gl_inside(gl, mask, CI(32, 3))
    @test length(gl) == 136
    @test length(expected) == 88
    @test length(WWFS.walk_gl_inside(gl, mask, CI(35, 40))) == 38
    for start in gl
        (; inside, outside) = WWFS.walk_gl(gl, mask, start)
        @test Set(inside) == Set(expected)
        @test allunique(inside)
        @test is_d4(inside)
        @test Set(outside) == Set(gl)
        @test allunique(outside)
        @test first(outside) == start
        unculled = WWFS.walk_gl(gl, mask, start; cull_inside=false)
        @test sort(length.(unculled.inside)) == [38, 88]
        @test all(allunique, unculled.inside)
        @test all(is_d4, unculled.inside)
        @test unculled.outside == outside
    end
end

@testset "Two inside GL walks sharing a junction" begin
    mask = falses(9, 9)
    mask[2:5, 2:5] .= true
    mask[5:7, 5:7] .= true
    gl = boundary(mask)
    expected = WWFS.walk_gl_inside(gl, mask, CI(2, 2))
    @test length(expected) == 12
    @test length(WWFS.walk_gl_inside(gl, mask, CI(7, 7))) == 8
    for start in gl
        (; inside, outside) = WWFS.walk_gl(gl, mask, start)
        @test Set(inside) == Set(expected)
        @test allunique(inside)
        @test is_d4([inside; first(inside)])
        @test Set(outside) == Set(gl)
        unculled = WWFS.walk_gl(gl, mask, start; cull_inside=false)
        @test sort(length.(unculled.inside)) == [8, 12]
        @test all(allunique, unculled.inside)
        @test all(path -> is_d4([path; first(path)]), unculled.inside)
        @test intersect(Set.(unculled.inside)...) == Set([CI(5, 5)])
        @test unculled.outside == outside
    end
end

@testset "Inside GL walk ties" begin
    mask = falses(9, 9)
    mask[2:4, 2:4] .= true
    mask[5:7, 5:7] .= true
    gl = boundary(mask)
    for start in (CI(2, 2), CI(7, 7))
        unculled = WWFS.walk_gl(gl, mask, start; cull_inside=false)
        @test length.(unculled.inside) == [8, 8]
        @test WWFS.walk_gl(gl, mask, start).inside == first(unculled.inside)
    end
end

@testset "Isolated GL cell and disconnected input" begin
    mask = falses(5, 5)
    mask[2, 2] = true
    gl = [CI(2, 2)]
    @test WWFS.walk_gl_outside(gl, mask) == gl
    @test WWFS.walk_gl(gl, mask) == (; inside=CI{2}[], outside=gl)
    @test WWFS.walk_gl(gl, mask; cull_inside=false) == (; inside=Vector{CI{2}}[], outside=gl)
    @test_throws ErrorException WWFS.walk_gl_inside(gl, mask)
    mask[4, 4] = true
    @test_throws ErrorException WWFS.walk_gl_outside([gl; CI(4, 4)], mask)
end

@testset "Inside GL halves stop before endpoint U-turns" begin
    mask = falses(7, 7)
    mask[3:5, 2:5] .= true
    mask[2, 5] = true  # Another outlet's GL allows the old group-local trace to turn back.
    gl = [CI(3, j) for j in 2:5]
    half = WWFS.trace_inside(Set(gl), mask, CI(3, 3), :right; previous=CI(3, 2))
    @test half == [CI(3, j) for j in 3:5]
    @test allunique(half)
    for start in gl
        path = WWFS.walk_gl_inside(gl, mask, start)
        @test Set(path) == Set(gl)
        @test allunique(path)
        @test is_d4(path)
    end
end

@testset "Classify interior using the complete GL boundary" begin
    mask = falses(8, 8)
    mask[2:5, 2:5] .= true
    mask[6, 6] = true
    gl = [CI(6, 6)]
    full_boundary = Set(boundary(mask))
    before = copy(mask), copy(gl), copy(full_boundary)
    @test WWFS.walk_gl(gl, mask).inside == gl  # Group-local classification is insufficient.
    @test !WWFS.in_region(CI(5, 5), full_boundary, mask, :inside)
    result = WWFS.walk_gl(gl, mask; boundary=full_boundary)
    @test isempty(result.inside)
    @test result.outside == gl
    @test isempty(WWFS.walk_gl(gl, mask; boundary=full_boundary, cull_inside=false).inside)
    @test (mask, gl, full_boundary) == before
end
end
