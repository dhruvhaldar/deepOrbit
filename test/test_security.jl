using Test
using DeepOrbit

@testset "Security Tests" begin
    lats = [0.0, 10.0]
    lons = [0.0, 10.0]

    @testset "Filename Extension Validation" begin
        # Valid extension
        @test_nowarn generate_ground_track_svg(lats, lons, "test_output.svg")
        @test_nowarn generate_ground_track_svg(lats, lons, "test_output.SVG") # Case insensitive

        # Invalid extensions
        @test_throws ArgumentError generate_ground_track_svg(lats, lons, "test_output.txt")
        @test_throws ArgumentError generate_ground_track_svg(lats, lons, "test_output")
        @test_throws ArgumentError generate_ground_track_svg(lats, lons, "test_output.svg.txt")

        # Cleanup
        rm("test_output.svg", force=true)
        rm("test_output.SVG", force=true)
    end
end
