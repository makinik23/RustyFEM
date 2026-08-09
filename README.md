# RustyFEM

RustyFEM is a finite element method solver written in Rust. The
current implementation focuses on static 2D analysis and keeps the numerical
steps visible: element matrices, global assembly, constraints, loads, solving,
and postprocessing.

## Documentation

RustyFEM currently implements four 2D finite elements:

- 2D truss
- 2D Euler-Bernoulli beam
- Plane-stress T3 triangle
- Plane-stress Q4 quadrilateral

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

The dense LU solver is used by default. Use the sparse CSR solver with
preconditioned Conjugate Gradient:

```bash
cargo run -- --space 2d --solver sparse \
  --cg-tolerance 1e-10 \
  --cg-max-iterations 1000
```

The program reads materials, sections, nodes, constraints, elements, and loads.
Loads can be entered as nodal loads, uniform line loads on beam elements, edge
tractions on plane-stress elements, body forces on plane-stress elements, or
self-weight on plane-stress elements.
Enter `done` at the end of each group.

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

The beam section and element input formats are:

```text
section> beam ID AREA I HEIGHT
element> beam ID NODE_1 NODE_2 MATERIAL_ID SECTION_ID
```

During the session, the program asks you to provide the model in several
groups. You must enter:

- material properties: Young's modulus, Poisson's ratio, and density;
- nodes: node ID and `x`, `y` coordinates;
- constraints: node ID, degree of freedom (`Ux`, `Uy`, or `Rz`), and prescribed displacement;
- sections: section type, section ID, and section properties;
- elements: element type, connectivity, material ID, and section ID;
- loads: nodal loads or supported element loads.

Enter `done` after each group to continue to the next prompt. For this beam,
the interactive session looks like this:

```text
material> 1 200 0.3 1
material> done
section> beam 2 1 2 2
section> done
node> 1 0 0
node> 2 1 0
node> done
constraint> 1 Ux 0
constraint> 1 Uy 0
constraint> 1 Rz 0
constraint> done
element> beam 10 1 2 1 2
element> done
load> 2 Uy -12
load> done
```

Plane-stress Q4 quadrilaterals can be entered with nodes ordered
counterclockwise around the element:

```text
q4 ID NODE_1 NODE_2 NODE_3 NODE_4 MATERIAL_ID SECTION_ID
```

The natural node order is bottom-left, bottom-right, top-right, top-left in
the reference square. The solver rejects inverted or degenerate Q4 elements by
checking the Jacobian determinant at the 2x2 Gauss points.

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

## Uniform Beam Line Loads

Uniform line loads can be applied to beam elements with:

```text
beam_uniform ELEMENT_ID local|global QX QY
```

`QX` and `QY` are force per unit beam length. With `local`, `QX` acts along the
beam axis and `QY` acts in the local transverse direction. With `global`, the
components are interpreted in the model x/y axes and transformed to the beam's
local axes during load-vector assembly.

For the beam example above, a downward uniform local load of `-12 N/m` is:

```text
load> beam_uniform 10 local 0 -12
```

The solver converts this load into the consistent equivalent nodal load vector.
Beam end-force and section-response recovery also include the fixed-end force
contribution from the uniform load.

## Plane-Stress Edge Tractions

Uniform tractions can be applied to one edge of a T3 or Q4 plane-stress element with:

```text
edge_traction ELEMENT_ID NODE_A NODE_B local|global TX TY
```

`NODE_A` and `NODE_B` must be the two nodes of one boundary edge of the
selected element. `TX` and `TY` are force per unit area. The solver multiplies
them by the plane-stress thickness and edge length, then distributes half of
the resulting force to each edge node.

With `global`, the components are interpreted in the model x/y axes. With
`local`, `TX` acts from `NODE_A` toward `NODE_B`, and `TY` acts in the edge
normal direction.

For example, this applies a downward global traction on edge `2-3` of element
`20`:

```text
load> edge_traction 20 2 3 global 0 -1000
```

## Plane-Stress Body Forces

Uniform body forces can be applied over the full area and thickness of a T3 or
Q4 plane-stress element with:

```text
body_force ELEMENT_ID global BX BY
```

`BX` and `BY` are force per unit volume in the model x/y axes. The solver
multiplies them by the element area and plane-stress thickness, then converts
the result to consistent nodal loads. For T3 this is one third of the total
force per node; for rectangular Q4 elements this is one quarter per node.

For example, this applies a downward body force to plane-stress element `20`:

```text
load> body_force 20 global 0 -9810
```

This is an explicitly supplied body-force intensity. Use `self_weight` when the
load should be calculated from material density.

## Plane-Stress Self-Weight Loads

Self-weight loads can be applied over the full area and thickness of a T3 or Q4
plane-stress element with:

```text
self_weight ELEMENT_ID global AX AY
```

`AX` and `AY` are acceleration components in the model x/y axes. The solver
multiplies them by the element material density, area, and plane-stress
thickness, then converts the result to consistent nodal loads. For T3 this is
one third of the total force per node; for rectangular Q4 elements this is one
quarter per node.

For example, this applies gravity in the negative global y direction:

```text
load> self_weight 20 global 0 -9.81
```

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

### T3 vs Q4 on the same grid

The comparison below uses the same rectangular grid nodes, supports, material,
thickness, and nodal right-edge loads for both element families. T3 splits each
rectangular cell into two triangles; Q4 uses one bilinear quadrilateral per
cell. The measured displacement is still the vertical displacement of the
middle node on the right edge.

| Mesh | Nodes | T3 elements | Q4 elements | T3 displacement [m] | T3 relative error | Q4 displacement [m] | Q4 relative error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `4x2` | `15` | `16` | `8` | `-0.473706556570` | `88.1573%` | `-1.169360954824` | `70.7660%` |
| `8x4` | `45` | `64` | `32` | `-1.371713039702` | `65.7072%` | `-2.491608327189` | `37.7098%` |
| `16x8` | `153` | `256` | `128` | `-2.700297683096` | `32.4926%` | `-3.482652929698` | `12.9337%` |
| `32x16` | `561` | `1024` | `512` | `-3.579778463598` | `10.5055%` | `-3.871268896729` | `3.2183%` |
| `64x32` | `2145` | `4096` | `2048` | `-3.901091199391` | `2.4727%` | `-3.983901623090` | `0.4025%` |

The same benchmark also compares the maximum recovered von Mises stress over
all plane-stress elements. T3 stress is constant inside each element. Q4 stress
is reported with two recovery modes: `GAUSS`, which uses the four
`2x2` integration points directly, and `CORNER/BILIN`, which bilinearly
extrapolates those Gauss-point stresses to the four element corners.

The Euler-Bernoulli stress reference uses the maximum bending moment at the
fixed end:

```text
M_max = |P| L = 10 N m
sigma_max = M_max (H / 2) / I = 60 Pa
von_mises_max = |sigma_max| = 60 Pa
```

At the outer fiber where bending stress is maximal, the rectangular-section
shear stress is zero, so the plane-stress von Mises value is equal to the
absolute bending stress.

| Mesh | Analytical max VM [Pa] | T3 max VM [Pa] | T3 error | Q4 GAUSS max VM [Pa] | Q4 GAUSS error | Q4 CORNER max VM [Pa] | Q4 CORNER error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `4x2` | `60.000000000000` | `13.027176757921` | `78.2880%` | `20.145515222150` | `66.4241%` | `30.959768052130` | `48.4004%` |
| `8x4` | `60.000000000000` | `25.515325713544` | `57.4745%` | `35.328832071679` | `41.1186%` | `45.928419651943` | `23.4526%` |
| `16x8` | `60.000000000000` | `41.835587401302` | `30.2740%` | `48.231110614778` | `19.6148%` | `54.680152786864` | `8.8664%` |
| `32x16` | `60.000000000000` | `53.073556709693` | `11.5441%` | `55.351625518350` | `7.7473%` | `59.095238081485` | `1.5079%` |
| `64x32` | `60.000000000000` | `59.512381474695` | `0.8127%` | `60.481008054251` | `0.8017%` | `63.322481920185` | `5.5375%` |

The regular comparison test covers the first four rows for both displacement
and stress:

```bash
cargo test t3_and_q4_cantilever_response_comparison_across_meshes -- --nocapture
```

The `64x32` stress row is generated with an ignored sparse benchmark:

```bash
cargo test t3_and_q4_cantilever_64x32_sparse_response_comparison -- --ignored --nocapture
```

The dense `64x32` comparison is also available, but it is much slower:

```bash
cargo test t3_and_q4_cantilever_64x32_dense_mesh_comparison -- --ignored --nocapture
```

## Dense vs Sparse Solver Benchmark

The same rectangular T3 cantilever models were solved with both the dense LU
solver and the sparse CSR solver using preconditioned Conjugate Gradient with
the Jacobi preconditioner. The measurements include stiffness assembly,
boundary-condition application, and the complete linear solve.

Benchmark command:

```bash
cargo bench --bench solver_benchmark -- \
  --sample-size 100 \
  --measurement-time 5 \
  --warm-up-time 5
```

Results were collected with Criterion in release mode. The values below are
the median time reported by Criterion on the benchmark machine.

| Mesh | Dense LU | Sparse PCG + Jacobi | Relative speed |
| --- | ---: | ---: | ---: |
| `4x2` | `11.777 us` | `25.990 us` | `0.45x` |
| `8x4` | `85.919 us` | `147.440 us` | `0.58x` |
| `16x8` | `2.1291 ms` | `977.280 us` | `2.18x` |
| `32x16` | `90.172 ms` | `8.0671 ms` | `11.18x` |
| `64x32` | `6.5227 s` | `79.985 ms` | `81.55x` |

The sparse solver is slower for small systems because COO/CSR assembly and
iterative-solver setup add fixed overhead. Starting at `16x8`, the reduced
storage and sparse matrix-vector products outweigh that overhead. For the
`64x32` mesh, the sparse solver is approximately 82 times faster than dense
LU. The largest dense measurements exceeded Criterion's five-second target
budget, but all configured samples were collected.

## Development Checks

Run the standard checks before committing changes:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```
