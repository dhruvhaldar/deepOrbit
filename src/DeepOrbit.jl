module DeepOrbit

using LinearAlgebra

# Constants
const Re = 6378.137          # Earth's equatorial radius [km]
const mu = 398600.4418       # Earth's gravitational parameter [km^3/s^2]
const J2 = 0.00108262668     # Second zonal harmonic
const omega_e = 7.2921159e-5 # Earth's rotation rate [rad/s]

# Precomputed constants
const J2_coef = 1.5 * J2 * mu * Re^2

export Re, mu, J2, omega_e, propagate, eci_to_geodetic, generate_ground_track_svg

include("Visualizations.jl")

"""
    state_derivative(t, state)

Computes the derivative of the state vector (position and velocity).
State is [x, y, z, vx, vy, vz].
"""
function state_derivative(t, state)
    r_vec = state[1:3]
    v_vec = state[4:6]
    
    r_norm = norm(r_vec)
    r2 = r_norm^2
    r3 = r_norm^3
    
    x, y, z = r_vec
    
    # Precompute common terms
    c1 = -mu / r3
    c2 = J2_coef / (r2 * r3)
    
    z2_r2 = (z / r_norm)^2
    
    ax = c1 * x + c2 * x * (5 * z2_r2 - 1)
    ay = c1 * y + c2 * y * (5 * z2_r2 - 1)
    az = c1 * z + c2 * z * (5 * z2_r2 - 3)
    
    return [v_vec[1], v_vec[2], v_vec[3], ax, ay, az]
end

"""
    state_derivative!(dstate, t, state)

In-place version of state_derivative. Updates dstate.
"""
function state_derivative!(dstate, t, state)
    @inbounds begin
        x = state[1]
        y = state[2]
        z = state[3]
        vx = state[4]
        vy = state[5]
        vz = state[6]
    end

    r2 = x*x + y*y + z*z

    # Optimization: Precompute inverse square distance to replace divisions with multiplications
    inv_r2 = 1.0 / r2
    r_norm = sqrt(r2)
    inv_r3 = inv_r2 / r_norm

    # J2 Perturbation terms
    c1 = -mu * inv_r3
    c2 = J2_coef * inv_r3 * inv_r2

    z2_r2 = (z * z) * inv_r2
    term_z = 5 * z2_r2

    # Factor out common terms
    # Optimization: reused common term for ax and ay, and simplified az
    common_term = c1 + c2 * (term_z - 1)
    ax = x * common_term
    ay = y * common_term
    az = z * (common_term - 2 * c2)

    @inbounds begin
        dstate[1] = vx
        dstate[2] = vy
        dstate[3] = vz
        dstate[4] = ax
        dstate[5] = ay
        dstate[6] = az
    end
    return nothing
end

"""
    rk4_step(f, t, y, dt)

Performs a single Runge-Kutta 4 step.
"""
function rk4_step(f, t, y, dt)
    k1 = f(t, y)
    k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = f(t + dt, y + dt * k3)
    
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
end

"""
    rk4_step!(next_state, f!, t, y, dt, k1, k2, k3, k4, temp_state)

In-place Runge-Kutta 4 step using pre-allocated buffers.
Optimized with explicit loops to avoid broadcast overhead on small vectors.
"""
function rk4_step!(next_state, f!, t, y, dt, k1, k2, k3, k4, temp_state)
    # k1 = f(t, y)
    f!(k1, t, y)

    dt_half = 0.5 * dt
    dt_sixth = dt / 6.0

    # k2 = f(t + 0.5*dt, y + 0.5*dt*k1)
    # Optimization: @simd allows vectorization of the loop
    @inbounds @simd for i in eachindex(y)
        temp_state[i] = y[i] + dt_half * k1[i]
    end
    f!(k2, t + dt_half, temp_state)

    # k3 = f(t + 0.5*dt, y + 0.5*dt*k2)
    @inbounds @simd for i in eachindex(y)
        temp_state[i] = y[i] + dt_half * k2[i]
    end
    f!(k3, t + dt_half, temp_state)

    # k4 = f(t + dt, y + dt*k3)
    @inbounds @simd for i in eachindex(y)
        temp_state[i] = y[i] + dt * k3[i]
    end
    f!(k4, t + dt, temp_state)

    # y_next = y + (dt/6)*(k1 + 2k2 + 2k3 + k4)
    @inbounds @simd for i in eachindex(y)
        next_state[i] = y[i] + dt_sixth * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
    end
    return nothing
end

"""
    propagate(initial_state, t_span, dt)

Propagates the orbit from `t_span[1]` to `t_span[2]` with time step `dt`.
Returns a tuple (times, states).
"""
function propagate(initial_state, t_span, dt)
    t0, tf = t_span

    # Security Validation
    if length(initial_state) != 6
        throw(ArgumentError("Initial state must have exactly 6 elements (pos + vel). Got $(length(initial_state))"))
    end

    if dt <= 0
        throw(ArgumentError("Time step dt must be positive. Got $dt"))
    end

    if tf < t0
        throw(ArgumentError("Time span must be increasing (tf >= t0). Got t0=$t0, tf=$tf"))
    end

    estimated_steps = (tf - t0) / dt
    MAX_STEPS = 10_000_000

    if estimated_steps > MAX_STEPS
        throw(ArgumentError("Simulation requires too many steps ($(round(estimated_steps)). > $MAX_STEPS). Increase dt or decrease t_span to prevent memory exhaustion."))
    end

    times = t0:dt:tf
    n_steps = length(times)
    
    # Pre-allocate states array: Allocate all vectors at once to avoid allocation in the loop
    states = [Vector{Float64}(undef, 6) for _ in 1:n_steps]
    copyto!(states[1], initial_state)
    
    # Reusable buffers for RK4 calculations
    k1 = Vector{Float64}(undef, 6)
    k2 = Vector{Float64}(undef, 6)
    k3 = Vector{Float64}(undef, 6)
    k4 = Vector{Float64}(undef, 6)
    temp_state = Vector{Float64}(undef, 6)
    
    @inbounds for i in 1:(n_steps - 1)
        t = times[i]
        # Direct write to next state vector in the array, using current state from array
        rk4_step!(states[i+1], state_derivative!, t, states[i], dt, k1, k2, k3, k4, temp_state)
    end
    
    return times, states
end

"""
    eci_to_geodetic(state, t)

Converts ECI position to Geodetic coordinates (Latitude, Longitude, Altitude).
Assumes t=0 corresponds to GMST=0 (Greenwich Meridian aligned with X-axis).
"""
function eci_to_geodetic(state, t)
    x = state[1]
    y = state[2]
    z = state[3]
    # OPTIMIZATION: Manually compute norm to avoid allocating [x, y, z] vector
    r = sqrt(x*x + y*y + z*z)
    
    # Calculate Right Ascension and Declination
    declination = asin(z / r)
    right_ascension = atan(y, x)
    
    # Calculate Greenwich Mean Sidereal Time (GMST) approximation
    # theta_gmst = theta_gmst0 + omega_e * t
    # Assuming theta_gmst0 = 0 at t=0
    theta_gmst = omega_e * t
    
    # Longitude = Right Ascension - GMST
    longitude = right_ascension - theta_gmst
    
    # Normalize longitude to [-pi, pi]
    longitude = (longitude + pi) % (2 * pi) - pi
    
    # Latitude (Geocentric) is approximately equal to declination for spherical earth
    # For geodetic latitude on ellipsoid, it's more complex, but for LEO viz, spherical is a good start.
    # However, let's just stick to Geocentric latitude (declination) for now as "Latitude".
    latitude = declination
    
    altitude = r - Re
    
    return rad2deg(latitude), rad2deg(longitude), altitude
end

end # module DeepOrbit
