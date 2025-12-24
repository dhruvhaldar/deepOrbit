module DeepOrbit

using LinearAlgebra

# Constants
const Re = 6378.137          # Earth's equatorial radius [km]
const mu = 398600.4418       # Earth's gravitational parameter [km^3/s^2]
const J2 = 0.00108262668     # Second zonal harmonic
const omega_e = 7.2921159e-5 # Earth's rotation rate [rad/s]

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
    c2 = 1.5 * J2 * mu * Re^2 / (r2 * r3)
    
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
    x = state[1]
    y = state[2]
    z = state[3]
    vx = state[4]
    vy = state[5]
    vz = state[6]

    r2 = x*x + y*y + z*z
    r_norm = sqrt(r2)
    r3 = r2 * r_norm

    # J2 Perturbation terms
    c1 = -mu / r3
    c2 = 1.5 * J2 * mu * Re^2 / (r2 * r3)

    z2_r2 = (z / r_norm)^2
    term_z = 5 * z2_r2

    ax = c1 * x + c2 * x * (term_z - 1)
    ay = c1 * y + c2 * y * (term_z - 1)
    az = c1 * z + c2 * z * (term_z - 3)

    dstate[1] = vx
    dstate[2] = vy
    dstate[3] = vz
    dstate[4] = ax
    dstate[5] = ay
    dstate[6] = az
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
"""
function rk4_step!(next_state, f!, t, y, dt, k1, k2, k3, k4, temp_state)
    # k1 = f(t, y)
    f!(k1, t, y)

    # k2 = f(t + 0.5*dt, y + 0.5*dt*k1)
    @. temp_state = y + 0.5 * dt * k1
    f!(k2, t + 0.5 * dt, temp_state)

    # k3 = f(t + 0.5*dt, y + 0.5*dt*k2)
    @. temp_state = y + 0.5 * dt * k2
    f!(k3, t + 0.5 * dt, temp_state)

    # k4 = f(t + dt, y + dt*k3)
    @. temp_state = y + dt * k3
    f!(k4, t + dt, temp_state)

    # y_next = y + (dt/6)*(k1 + 2k2 + 2k3 + k4)
    @. next_state = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
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

    times = collect(t0:dt:tf)
    n_steps = length(times)
    
    # Pre-allocate states array
    states = Vector{Vector{Float64}}(undef, n_steps)
    states[1] = copy(initial_state)
    
    # Reusable buffers to avoid allocation in loop
    current_state = copy(initial_state)
    next_state = Vector{Float64}(undef, 6)
    k1 = Vector{Float64}(undef, 6)
    k2 = Vector{Float64}(undef, 6)
    k3 = Vector{Float64}(undef, 6)
    k4 = Vector{Float64}(undef, 6)
    temp_state = Vector{Float64}(undef, 6)
    
    for i in 1:(n_steps - 1)
        t = times[i]
        rk4_step!(next_state, state_derivative!, t, current_state, dt, k1, k2, k3, k4, temp_state)

        # Save result (must allocate a new vector for the results array)
        states[i+1] = copy(next_state)

        # Update current_state for next step (copy without allocation)
        copyto!(current_state, next_state)
    end
    
    return times, states
end

"""
    eci_to_geodetic(state, t)

Converts ECI position to Geodetic coordinates (Latitude, Longitude, Altitude).
Assumes t=0 corresponds to GMST=0 (Greenwich Meridian aligned with X-axis).
"""
function eci_to_geodetic(state, t)
    x, y, z = state[1:3]
    r = norm([x, y, z])
    
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
