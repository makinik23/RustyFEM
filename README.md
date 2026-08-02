# RustyFEM

RustyFEM is a finite element method solver written in Rust. The
current implementation focuses on static 2D analysis and keeps the numerical
steps visible: element matrices, global assembly, constraints, loads, solving,
and postprocessing.

## Documentation

RustyFEM currently implements three 2D finite elements:

- 2D truss
- 2D Euler-Bernoulli beam
- Plane-stress T3 triangle

The [engineering assumptions](docs/engineering-assumptions.md) document
explains the material model, static analysis assumptions, degrees of freedom,
and limitations of each formulation. Continue there to understand what the
solver represents mechanically before reading the implementation.

The examples below show the theory in action: a beam is compared with analytical formulas for its deflection and bending moment,
while a sequence of progressively finer triangular meshes demonstrates how the numerical solution converges toward an analytical solution.

## Requirements

- Rust 1.89 or newer
- Cargo

## Run The Interactive Solver

Start a 2D interactive session with:

```bash
cargo run -- --space 2d
```

The program reads one material, then nodes, constraints, elements, and nodal
loads. Enter `done` at the end of each group.

## Beam Example

The following example is a cantilever beam with:

- length `L = 1 m`
- Young's modulus `E = 200 Pa`
- density `rho = 1 kg/m^3`
- cross-sectional area `A = 1 m^2`
- second moment of area `I = 2 m^4`
- section height `h = 2 m`
- downward tip load `P = -12 N`
- the first node fully fixed

The beam input format is:

```text
beam ID NODE_1 NODE_2 AREA I HEIGHT
```

During the session, the program asks you to provide the model in several
groups. You must enter:

- material properties: Young's modulus, Poisson's ratio, and density;
- nodes: node ID and `x`, `y` coordinates;
- constraints: node ID, degree of freedom (`Ux`, `Uy`, or `Rz`), and prescribed displacement;
- elements: element type, connectivity, and section properties;
- loads: node ID, degree of freedom, and load value.

Enter `done` after each group to continue to the next prompt. For this beam,
the interactive session looks like this:

```text
material> 200 0.3 1
node> 1 0 0
node> 2 1 0
node> done
constraint> 1 Ux 0
constraint> 1 Uy 0
constraint> 1 Rz 0
constraint> done
element> beam 10 1 2 1 2 2
element> done
load> 2 Uy -12
load> done
```

For a cantilever with a point load at the free end, the analytical solution is:

```text
v(L)     = P L^3 / (3 E I) = -0.010 m
theta(L) = P L^2 / (2 E I) = -0.015 rad
```

The numerical beam result agrees with the analytical values to within machine precision (floating-point round-off).

The recovered bending moment and fiber stress also agree with the analytical
relations:

```text
M(x)       = P (L - x)       N m
kappa(x)  = M(x) / (E I)     1/m
sigma(y)  = N / A - M(x) y / I Pa
```

At the fixed end, `M(0) = -12 N m` and `y = +/- h/2 = +/- 1 m`, so the
recovered fiber stresses are `+6 Pa` and `-6 Pa`.

## T3 Triangle Mesh Benchmark

The T3 benchmark models a rectangular plane-stress cantilever:

- rectangle dimensions: `L = 10 m`, `H = 1 m`
- thickness: `t = 1 m`
- material: `E = 1000 Pa`, `nu = 0.3 [-]`
- left edge: `Ux = Uy = 0`
- total downward load on the right edge: `P = -1 N`
- every rectangular cell is split into two T3 triangles

For comparison, the Euler-Bernoulli beam reference uses:

```text
I = t H^3 / 12 = 1/12 m^4
delta_analytical = P L^3 / (3 E I) = -4.0 m
```

The measured displacement is the vertical displacement of the middle node on
the right edge and is reported in `m`. The table shows convergence as the
triangular mesh is refined. Relative error is shown as a
percentage.
The dimensions in the first column are the number of rectangular cells in the
`x` and `y` directions; each cell contains two T3 elements.

| Mesh | Numerical displacement [m] | Analytical displacement [m] | Relative error |
| --- | ---: | ---: | ---: |
| `4x2` | `-0.473706556570` | `-4.000000000000` | `88.1573%` |
| `8x4` | `-1.371713039702` | `-4.000000000000` | `65.7072%` |
| `16x8` | `-2.700297683096` | `-4.000000000000` | `32.4926%` |
| `32x16` | `-3.579778463598` | `-4.000000000000` | `10.5055%` |
| `64x32` | `-3.901091199391` | `-4.000000000000` | `2.4727%` |

## Development Checks

Run the standard checks before committing changes:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```
