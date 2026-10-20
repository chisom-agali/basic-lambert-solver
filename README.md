# lambert-solver

A basic 2D Lambert solver in Python. You give it two position vectors and a time of flight, and it gives back the departure and arrival velocities for the orbit that connects them. I kept this version simple on purpose so I could get the fundamentals down first. It only does one short-way, prograde transfer with zero full revolutions.

## Setup

```
pip install -r requirements.txt
```

## Run

```
python interplanetary_transfer_basic.py
```

This solves one example (Earth orbit, r1 = [7000, 0] km to r2 = [0, 8000] km in 1 hour), prints the velocities, and saves a plot of the geometry as `lambert_geometry.png`.

## Units

Distance is in km, time in seconds, velocity in km/s, and mu in km^3/s^2.

## How it works

**Stumpff functions.** `stumpff_C` and `stumpff_S` are the pieces of the universal variable method. They let the same equations work whether the trial orbit is an ellipse or a hyperbola, so I don't need separate cases.

**The solver.** First it finds the angle between r1 and r2 and builds a geometry term A. Then it guesses a value of z (which controls the shape of the trial orbit) and checks how far the time of flight is from what I asked for. It scans a range of z values to find where that error changes sign, then uses bisection to zero in on the answer.

**Getting the velocities.** Once z is solved, I use the Lagrange f and g coefficients to turn the geometry into v1 and v2.

**The example.** `main()` runs the Earth orbit case and plots r1, r2, and the two velocity vectors. The velocity arrows are scaled up so they're visible, so their length isn't to scale with the positions.

## Check

I propagated the departure velocity forward for 1 hour with a numerical integrator and it landed on r2 to within a few millimeters. Energy and angular momentum also match at both ends, which they should since it's a single Keplerian orbit.

## Limits

The transfer angle comes from `acos`, so it's always between 0 and 180 degrees. That means it can't handle prograde transfers that go more than halfway around, and it fails if r1 and r2 are lined up (0 or 180 degrees). Multi-revolution solutions, retrograde, and 3D are things I'd add later.
