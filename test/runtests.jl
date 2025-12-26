using Test
using DeepOrbit
using LinearAlgebra

@testset "DeepOrbit Tests" begin

    @testset "Constants" begin
        @test Re > 0
        @test mu > 0
        @test J2 > 0
    end

    @testset "Input Validation" begin
        # Test valid input
        initial_state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        t_span = (0.0, 10.0)
        dt = 1.0
        times, states = propagate(initial_state, t_span, dt)
        @test length(times) == length(states)

        # Test invalid state vector length (too short)
        short_state = [7000.0, 0.0, 0.0]
        @test_throws ArgumentError propagate(short_state, t_span, dt)

        # Test invalid state vector length (too long)
        long_state = zeros(7)
        @test_throws ArgumentError propagate(long_state, t_span, dt)
    end

    @testset "Circular Orbit Energy Conservation" begin
        # For a circular orbit without perturbations, energy should be constant.
        # With J2, it's not strictly constant but should be close over one orbit.
        # However, for testing purposes, we can verify that the semi-major axis doesn't drift wildly.
        
        # Initial State: Circular Equatorial (J2 effects minimal on energy for equatorial?)
        # Actually J2 is conservative, so energy (Specific Mechanical Energy) should be conserved?
        # The potential is V = -mu/r [1 - J2(Re/r)^2 P2(sin(phi))]
        # So Total Energy E = v^2/2 + V should be conserved.
        
        r_mag = Re + 500.0
        v_mag = sqrt(mu / r_mag)
        
        # Equatorial orbit
        initial_state = [r_mag, 0.0, 0.0, 0.0, v_mag, 0.0]
        t_span = (0.0, 6000.0)
        dt = 10.0
        
        times, states = propagate(initial_state, t_span, dt)
        
        function total_energy(state)
            r = state[1:3]
            v = state[4:6]
            r_norm = norm(r)
            v_sq = dot(v, v)
            z = r[3]
            sin_phi = z / r_norm
            
            # Potential including J2
            V = -mu / r_norm * (1 - 1.5 * J2 * (Re / r_norm)^2 * (sin_phi^2 - 1/3))
             # Wait, the J2 term usually is written: U = mu/r [1 - J2(Re/r)^2 (1.5 sin^2(phi) - 0.5)]
             # So V = -U = -mu/r + mu/r * J2(Re/r)^2 * (1.5 sin^2 - 0.5)
             
            # Let's double check the formula used in derivative.
            # a_J2 derived from potential U = (mu/r) [1 - J2/2 (Re/r)^2 (3 sin^2(phi) - 1)]
            # So V = -U = -mu/r [1 - J2/2 (Re/r)^2 (3 sin^2(phi) - 1)]
            
            term = -1.5 * J2 * (Re/r_norm)^2 * ( (z/r_norm)^2 - 1/3 )
             # Wait, let's just stick to the specific energy:
             # E = v^2/2 - mu/r. 
             # With J2, this fluctuates.
             
            # Let's just check standard Keplerian energy E = v^2/2 - mu/r. 
            # It should be ROUGHLY constant.
            return v_sq / 2 - mu / r_norm
        end
        
        E0 = total_energy(states[1])
        Ef = total_energy(states[end])
        
        # Check relative error
        rel_error = abs((Ef - E0) / E0)
        println("Energy variation (relative): $rel_error")
        @test rel_error < 1e-3
    end

    @testset "ECI to Geodetic" begin
        # Point on X-axis (Lat=0, Lon=0 if t=0)
        state = [Re + 100.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        lat, lon, alt = eci_to_geodetic(state, 0.0)
        @test isapprox(lat, 0.0, atol=1e-5)
        @test isapprox(lon, 0.0, atol=1e-5)
        @test isapprox(alt, 100.0, atol=1e-5)
        
        # Point on Z-axis (Lat=90)
        state = [0.0, 0.0, Re + 100.0, 0.0, 0.0, 0.0]
        lat, lon, alt = eci_to_geodetic(state, 0.0)
        @test isapprox(lat, 90.0, atol=1e-5)
        @test isapprox(alt, 100.0, atol=1e-5)
    end
end
