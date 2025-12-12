module DeepOrbit

using LinearAlgebra

# Constants
const Re = 6378.137          # Earth's equatorial radius [km]
const mu = 398600.4418       # Earth's gravitational parameter [km^3/s^2]
const J2 = 0.00108262668     # Second zonal harmonic
const omega_e = 7.2921159e-5 # Earth's rotation rate [rad/s]

export Re, mu, J2, omega_e, propagate, eci_to_geodetic

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
    
    # J2 Perturbation terms
    # Formula from "Fundamentals of Astrodynamics and Applications", Vallado, or similar standard texts
    # a_J2_x = - (3/2) * J2 * (mu/r^2) * (Re/r)^2 * (x/r) * (1 - 5*(z/r)^2)
    # a_J2_y = - (3/2) * J2 * (mu/r^2) * (Re/r)^2 * (y/r) * (1 - 5*(z/r)^2)
    # a_J2_z = - (3/2) * J2 * (mu/r^2) * (Re/r)^2 * (z/r) * (3 - 5*(z/r)^2)
    
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
    propagate(initial_state, t_span, dt)

Propagates the orbit from `t_span[1]` to `t_span[2]` with time step `dt`.
Returns a tuple (times, states).
"""
function propagate(initial_state, t_span, dt)
    t0, tf = t_span
    times = collect(t0:dt:tf)
    n_steps = length(times)
    
    # Pre-allocate states array
    # Rows are steps, columns are state variables (6)
    states = Vector{Vector{Float64}}(undef, n_steps)
    states[1] = initial_state
    
    current_state = initial_state
    
    for i in 1:(n_steps - 1)
        t = times[i]
        current_state = rk4_step(state_derivative, t, current_state, dt)
        states[i+1] = current_state
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
