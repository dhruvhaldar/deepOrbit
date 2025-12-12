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

This will generate a `ground_track.csv` file containing the Time, Latitude, Longitude, and Altitude of the satellite, as well as a `ground_track.svg` file visualizing the orbit path.

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

## Technical Details

### Physics Model (J2 Perturbation)

The simulation uses a **Runge-Kutta 4th Order (RK4)** numerical integrator to propagate the satellite's state vector over time. The primary force acting on the satellite is Earth's gravity, modeled with the **J2 zonal harmonic perturbation**. This accounts for Earth's oblateness (bulge at the equator), which causes orbital orbital plane precession (nodal regression) and rotation of the perigee.

The acceleration vector `a` is computed as:

```julia
a = -mu/r^3 * r + a_J2
```

Where `a_J2` contains the perturbation terms derived from the specific potential energy.

### Constants & Variables

The following key constants are defined in `src/DeepOrbit.jl`:

- `Re` (**6378.137 km**): Earth's equatorial radius.
- `mu` (**398600.4418 km³/s²**): Earth's standard gravitational parameter.
- `J2` (**0.00108262668**): The second zonal harmonic coefficient representing Earth's oblateness.
- `omega_e` (**7.2921159e-5 rad/s**): Earth's mean rotation rate (used for converting Inertial -> Fixed coordinates).

### State Vector

The state vector is a 6-element array representing the satellite's kinematic state in the **Earth-Centered Inertial (ECI)** frame:

```julia
state = [x, y, z, vx, vy, vz]
```

- `x, y, z`: Position components in Kilometers (km).
- `vx, vy, vz`: Velocity components in Kilometers per second (km/s).

### Coordinate Systems

1.  **ECI (Earth-Centered Inertial):**

    - Origin: Earth's center of mass.
    - X-axis: Points towards the Vernal Equinox.
    - Z-axis: Points along the celestial noorth pole.
    - Used for: Numerical integration (physics).

2.  **ECEF (Earth-Centered Earth-Fixed) / Geodetic:**
    - Rotates with the Earth.
    - Used for: Ground track visualization (Latitude, Longitude).
    - The conversion in `eci_to_geodetic` accounts for Earth's rotation `omega_e * t` (Greenwich Mean Sidereal Time approximation).

### Visualization

The plotting logic (`src/Visualizations.jl`) generates an SVG file by performing an **Equirectangular projection**:

- **Longitude** (-180° to +180°) maps linearly to the X-axis.
- **Latitude** (-90° to +90°) maps linearly to the Y-axis.
- **Date Line Handling:** The plotter detects large jumps in longitude (> 300°) between consecutive points and breaks the SVG path (`M` move command) to prevent horizontal lines streaking across the map.

## Tests

To run the test suite:

```bash
julia --project=. test/runtests.jl
```
