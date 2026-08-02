# Engineering Assumptions

RustyFEM currently implements small-strain, linear static finite element
analysis in two dimensions. This document records the mechanical assumptions
behind the implemented formulations.

## Material Model

- The material is homogeneous, isotropic, and linearly elastic.
- Young's modulus `E` is finite and strictly positive.
- Poisson's ratio `nu` must satisfy `-1 < nu < 0.5`.
- Density is stored and validated, but it is not used by the current static
  formulation because inertia and body forces are not implemented.
- A `Model2D` currently contains one material shared by all its elements.
- Plasticity, damage, anisotropy, temperature dependence, and other nonlinear
  constitutive effects are outside the current scope.

## Analysis Assumptions

- The analysis is static and is based on the linear equilibrium equation
  `K u = f`.
- Displacements and rotations are assumed to be small.
- Geometric nonlinearity, large rotations, contact, buckling, dynamics, and
  time-dependent effects are not implemented.
- Loads are applied as nodal loads.
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

## Current Limitations

- The 3D analysis mode is recognized by the CLI but is not implemented yet.
- Body forces, inertia, dynamics, and non-nodal load distributions are not
  implemented.
