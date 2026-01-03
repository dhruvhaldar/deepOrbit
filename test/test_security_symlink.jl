using Test
using DeepOrbit

@testset "Security Tests - Symlinks" begin
    # Mock data
    lats = [0.0, 10.0]
    lons = [0.0, 10.0]

    # Create a target file
    target_file = "sensitive_target.txt"
    write(target_file, "Sensitive Content")

    # Create a symlink
    symlink_file = "symlink_test.svg"
    # Clean up if exists
    rm(symlink_file, force=true)

    symlink(target_file, symlink_file)

    @testset "Symlink Attack Prevention" begin
        # Verify it IS a symlink
        @test islink(symlink_file)

        # Verify function throws error
        @test_throws ArgumentError generate_ground_track_svg(lats, lons, symlink_file)

        # Verify content wasn't changed
        content = read(target_file, String)
        @test content == "Sensitive Content"
    end

    # Clean up
    rm(symlink_file, force=true)
    rm(target_file, force=true)
end
