using DeepOrbit
using DelimitedFiles

# Initial State (approx ISS)
# Altitude ~ 400 km
r_mag = Re + 400.0
v_mag = sqrt(mu / r_mag) # Circular orbit velocity

# Position on X-axis, Velocity on Y-axis (Equatorial orbit for simplicity, but let's incline it)
inclination = deg2rad(51.6) # ISS inclination
vx = 0.0
vy = v_mag * cos(inclination)
vz = v_mag * sin(inclination)

initial_state = [r_mag, 0.0, 0.0, vx, vy, vz]

# Simulation parameters
t_span = (0.0, 5400.0) # Approx 90 mins (one orbit)
dt = 10.0 # 10 seconds step

println("Starting simulation...")
times, states = propagate(initial_state, t_span, dt)
println("Simulation complete.")

# Process results
output_data = Matrix{Float64}(undef, length(times), 4) # Time, Lat, Lon, Alt

for i in 1:length(times)
    t = times[i]
    state = states[i]
    lat, lon, alt = eci_to_geodetic(state, t)
    
    output_data[i, 1] = t
    output_data[i, 2] = lat
    output_data[i, 3] = lon
    output_data[i, 4] = alt
end

# Save to CSV
output_file = "ground_track.csv"
writedlm(output_file, output_data, ',')
println("Results saved to $output_file")

# Print first and last few points
println("First 5 points:")
display(output_data[1:5, :])
println("\nLast 5 points:")
display(output_data[end-4:end, :])
