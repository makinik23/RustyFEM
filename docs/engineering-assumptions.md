# Engineering Assumptions

RustyFEM currently implements small-strain, linear static finite element
analysis in two dimensions. This document records the mechanical assumptions
behind the implemented formulations.

## Material Model

- The material is homogeneous, isotropic, and linearly elastic.
- Young's modulus `E` is finite and strictly positive.
- Poisson's ratio `nu` must satisfy `-1 < nu < 0.5`.
- Density is stored and validated. It is used by self-weight loads and is not
  otherwise used by the static stiffness formulation.
- A `Model2D` currently contains one material shared by all its elements.
- Plasticity, damage, anisotropy, temperature dependence, and other nonlinear
  constitutive effects are outside the current scope.

## Analysis Assumptions

- The analysis is static and is based on the linear equilibrium equation
  `K u = f`.
- Displacements and rotations are assumed to be small.
- Geometric nonlinearity, large rotations, contact, buckling, dynamics, and
  time-dependent effects are not implemented.
- Loads can be applied as nodal loads, uniform line loads on 2D beam elements,
  uniform edge tractions on T3/Q4 plane-stress elements, or uniform body
  forces and self-weight loads on T3/Q4 plane-stress elements.
- Uniform beam line loads are converted to consistent equivalent nodal loads
  before solving. They can be defined in the beam's local axis system or in the
  model's global x/y axis system.
- Plane-stress edge tractions are converted to equivalent nodal loads before
  solving. They can be defined in the global x/y axis system or in an
  edge-local system whose x-axis runs from the first specified edge node to the
  second.
- Plane-stress body forces are converted to equivalent nodal loads before
  solving. They are defined as force per unit volume in the global x/y axis
  system.
- Plane-stress self-weight loads are converted to equivalent nodal loads
  before solving. They are defined as acceleration in the global x/y axis
  system and are multiplied by the loaded element material density.
- Boundary conditions are prescribed nodal displacements, including zero
  displacement supports.
- The global stiffness matrix is assembled in the model's global degrees of
  freedom order.
- The current reference solver uses a dense direct LU solve.

## 2D Truss

- The truss is a two-node axial member with two translational DOFs per node:
  `Ux` and `Uy`.
- Linear Lagrange interpolation is used along the element axis.
- The element carries axial force and axial stress only.
- Bending, shear deformation, torsion, and nodal rotations are not represented.
- The cross-sectional area is constant along the element.

## 2D Euler-Bernoulli Beam

- The beam is a two-node Euler-Bernoulli frame element with `Ux`, `Uy`, and
  `Rz` at each node.
- Cubic Hermite interpolation is used for the transverse displacement and
  linear interpolation for the axial displacement.
- Plane cross-sections remain plane and normal to the neutral axis.
- Transverse shear deformation is neglected.
- The cross-sectional area and second moment of area are constant along the
  element.
- The formulation supports axial force, shear force, bending moment,
  curvature, and recovered top/bottom fiber stress.
- Uniform line loads on beams are included in equivalent nodal loads and in
  fixed-end force recovery for beam end forces and section moments.

## Plane-Stress T3 Triangle

- The T3 element is a three-node linear triangle with two translational DOFs
  per node: `Ux` and `Uy`.
- The formulation assumes plane stress, so the out-of-plane stress is zero:
  `sigma_z = 0`.
- The displacement field is linear inside the element.
- Strain and stress are constant within each T3 element.
- Thickness is constant and the element has no independent bending DOF.
- The element uses isotropic linear elasticity and can recover in-plane strain,
  stress, and von Mises stress.
- Uniform edge tractions can be applied to any one of the three T3 edges. The
  total edge force is `traction * thickness * edge_length` and is distributed
  equally to the two edge nodes.
- Uniform body forces can be applied over the whole T3 element. The total body
  force is `body_force * area * thickness` and is distributed equally to the
  three element nodes.
- Self-weight loads can be applied over the whole T3 element. The total force
  is `density * acceleration * area * thickness` and is distributed equally to
  the three element nodes.

## Plane-Stress Q4 Quadrilateral

- The Q4 element is a four-node bilinear quadrilateral with two translational
  DOFs per node: `Ux` and `Uy`.
- The formulation assumes plane stress, so the out-of-plane stress is zero:
  `sigma_z = 0`.
- The same bilinear shape functions interpolate geometry and displacement in
  the natural square `-1 <= xi <= 1`, `-1 <= eta <= 1`.
- The element stiffness is integrated with a full `2x2` Gauss rule.
- The Jacobian determinant must be finite and positive at all four Gauss
  points. This rejects inverted or degenerate node orderings.
- Strain and stress vary within the element. The common element-response API
  reports the Q4 response at the element center, `xi = 0`, `eta = 0`.
- Q4 stress recovery also supports Nastran-style output modes: `Center`,
  `Gauss` at the four `2x2` integration points, and `Corner`/`BILIN` by
  bilinear extrapolation from the Gauss points to the natural corners.
- Thickness is constant and the element has no independent bending DOF.
- Uniform edge tractions can be applied to any one of the four Q4 edges. The
  total edge force is `traction * thickness * edge_length` and is distributed
  equally to the two edge nodes.
- Uniform body forces and self-weight loads are integrated with the same `2x2`
  Gauss rule used for stiffness. For rectangular Q4 elements with uniform
  loading, the total force is distributed equally to the four element nodes.

## Current Limitations

- The 3D analysis mode is recognized by the CLI but is not implemented yet.
- Inertia, dynamics, nonuniform load distributions, and higher-order
  surface/edge tractions are not implemented.
