# DeepOrbit.jl

**DeepOrbit** is a high-performance Julia library for simulating satellite orbits in Low Earth Orbit (LEO). It provides accurate propagation of satellite state vectors accounting for the **J2 zonal harmonic perturbation**, which arises from the Earth's oblateness.

This project is designed for "deep tech" applications where precise orbital dynamics are required for analysis, visualization, or mission planning.

## Features

- **J2 Perturbation Physics:** Accurate modeling of the Earth's oblateness effect on satellite orbits.
- **RK4 Integration:** Robust Runge-Kutta 4th order numerical integrator for stable propagation.
- **Ground Track Calculation:** Built-in utilities to convert Earth-Centered Inertial (ECI) coordinates to Geodetic coordinates (Latitude, Longitude, Altitude).
- **High Performance:** Written in pure Julia for speed and efficiency.

## Installation

You can install the package by cloning this repository. Since the project relies only on Julia standard libraries, no additional package instantiation is required.

```bash
git clone https://github.com/your-username/DeepOrbit.jl.git
cd DeepOrbit
```

## Usage

### Running a Simulation

An example script is provided in `examples/simulate.jl` which simulates an ISS-like orbit (51.6° inclination) for one orbital period.

```bash
julia --project=. examples/simulate.jl
```

This will generate a `ground_track.csv` file containing the Time, Latitude, Longitude, and Altitude of the satellite.

### Library Usage

```julia
using DeepOrbit

# Define initial state (Position [km], Velocity [km/s])
r_mag = Re + 400.0
v_mag = sqrt(mu / r_mag)
initial_state = [r_mag, 0.0, 0.0, 0.0, v_mag, 0.0]

# Propagate for 1 hour (3600 seconds)
times, states = propagate(initial_state, (0.0, 3600.0), 10.0)

# Convert to Geodetic coordinates
lat, lon, alt = eci_to_geodetic(states[end], times[end])
println("Final Position: Lat=$lat, Lon=$lon, Alt=$alt")
```

## Tests

To run the test suite:

```bash
julia --project=. test/runtests.jl
```
