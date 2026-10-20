"""First we'll assume a simple 2-D Lambert solver, with zero full revolutions.

Units:
    distance: km
    time: seconds
    velocity: km/s
    mu: km^3/s^2

This version intentionally solves one short-way prograde transfer only.
More advanced concepts can be added
after the fundamentals are learned.
"""

import math
import numpy as np
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Stumpff functions
# -----------------------------------------------------------------------------
# These functions are part of the universal-variable Lambert method.
# They let the same equations work for elliptical and hyperbolic trial orbits.

def stumpff_C(z):
    """Return the C Stumpff function."""
    if abs(z) < 1e-6:
        return 0.5 - z / 24.0

    if z > 0:
        root_z = math.sqrt(z)
        return (1.0 - math.cos(root_z)) / z

    root_minus_z = math.sqrt(-z)
    return (math.cosh(root_minus_z) - 1.0) / (-z)


def stumpff_S(z):
    """Return the S Stumpff function."""
    if abs(z) < 1e-6:
        return 1.0 / 6.0 - z / 120.0

    if z > 0:
        root_z = math.sqrt(z)
        return (root_z - math.sin(root_z)) / root_z**3

    root_minus_z = math.sqrt(-z)
    return (math.sinh(root_minus_z) - root_minus_z) / root_minus_z**3


# -----------------------------------------------------------------------------
# Lambert solver
# -----------------------------------------------------------------------------

def solve_lambert_2d(r1, r2, time_of_flight, mu):
    """Solve one short-way, prograde 2-D Lambert problem.

    Parameters
    ----------
    r1 : array-like
        Starting position [x, y] in km.
    r2 : array-like
        Ending position [x, y] in km.
    time_of_flight : float
        Transfer time in seconds.
    mu : float
        Gravitational parameter in km^3/s^2.

    Returns
    -------
    v1, v2 : numpy arrays
        Departure and arrival velocities in km/s.
    """
    r1 = np.array(r1, dtype=float)
    r2 = np.array(r2, dtype=float)

    r1_size = np.linalg.norm(r1)
    r2_size = np.linalg.norm(r2)

    if r1_size == 0 or r2_size == 0:
        raise ValueError("Position vectors cannot be zero.")
    if time_of_flight <= 0:
        raise ValueError("Time of flight must be positive.")

    # The angle between the position vectors.
    cosine_angle = np.dot(r1, r2) / (r1_size * r2_size)
    cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
    transfer_angle = math.acos(cosine_angle)

    sine_angle = math.sin(transfer_angle)
    if abs(sine_angle) < 1e-8:
        raise ValueError("The position vectors are nearly collinear.")

    # A is a geometry term used in the universal-variable equations.
    A = sine_angle * math.sqrt(
        r1_size * r2_size / (1.0 - cosine_angle)
    )

    def time_error(z):
        """Difference between calculated and requested transfer time."""
        C = stumpff_C(z)
        S = stumpff_S(z)

        if C <= 0:
            return math.nan

        y = r1_size + r2_size + A * (z * S - 1.0) / math.sqrt(C)
        if y < 0:
            return math.nan

        x = math.sqrt(y / C)
        calculated_time = (x**3 * S + A * math.sqrt(y)) / math.sqrt(mu)
        return calculated_time - time_of_flight

    # Find a sign-changing interval for z.
    z_values = np.linspace(-4.0 * math.pi**2, 4.0 * math.pi**2, 2000)
    lower_z = None
    upper_z = None
    previous_z = None
    previous_error = None

    for z in z_values:
        error = time_error(z)

        if not math.isfinite(error):
            continue

        if previous_error is not None and error * previous_error < 0:
            lower_z = previous_z
            upper_z = z
            break

        previous_z = z
        previous_error = error

    if lower_z is None:
        raise RuntimeError("Could not find a Lambert solution bracket.")

    # Bisection: repeatedly cut the interval in half.
    lower_error = time_error(lower_z)

    for _ in range(100):
        z = 0.5 * (lower_z + upper_z)
        middle_error = time_error(z)

        if abs(middle_error) < 1e-6:
            break

        if lower_error * middle_error <= 0:
            upper_z = z
        else:
            lower_z = z
            lower_error = middle_error

    # Reconstruct the endpoint velocities using the f and g coefficients.
    C = stumpff_C(z)
    S = stumpff_S(z)
    y = r1_size + r2_size + A * (z * S - 1.0) / math.sqrt(C)

    f = 1.0 - y / r1_size
    g = A * math.sqrt(y / mu)
    g_dot = 1.0 - y / r2_size

    v1 = (r2 - f * r1) / g
    v2 = (g_dot * r2 - r1) / g

    return v1, v2


# -----------------------------------------------------------------------------
# Example problem
# -----------------------------------------------------------------------------

def main():
    """Run one example so we can inspect the solver output."""
    mu_earth = 398600.4418  # km^3/s^2

    r1 = np.array([7000.0, 0.0])
    r2 = np.array([0.0, 8000.0])
    time_of_flight = 3600.0  # seconds

    v1, v2 = solve_lambert_2d(r1, r2, time_of_flight, mu_earth)

    print("Starting position r1 (km):", r1)
    print("Ending position r2 (km):  ", r2)
    print("Departure velocity v1 (km/s):", v1)
    print("Arrival velocity v2 (km/s):  ", v2)

    # Plot the position and velocity vectors.
    plt.figure(figsize=(8, 8))

    plt.plot([0, r1[0]], [0, r1[1]], "bo-", label="Position r1")
    plt.plot([0, r2[0]], [0, r2[1]], "go-", label="Position r2")

    # The velocity vectors are scaled only so they are visible.
    velocity_scale = 500.0

    plt.quiver(
        r1[0], r1[1],
        v1[0] * velocity_scale,
        v1[1] * velocity_scale,
        color="red",
        angles="xy",
        scale_units="xy",
        scale=1,
        label="Departure velocity v1",
    )

    plt.quiver(
        r2[0], r2[1],
        v2[0] * velocity_scale,
        v2[1] * velocity_scale,
        color="purple",
        angles="xy",
        scale_units="xy",
        scale=1,
        label="Arrival velocity v2",
    )

    plt.xlabel("x position (km)")
    plt.ylabel("y position (km)")
    plt.title("Basic 2D Lambert Geometry")
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.savefig("lambert_geometry.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
