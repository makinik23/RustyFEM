# RustyFEM

RustyFEM is an educational finite element method solver written in Rust.

The project starts with a technically correct but deliberately small milestone:
a linear static solver for one-dimensional axial bar elements.

## Units and Indexing

All values use SI units:

- meters for length and displacement,
- square meters for area,
- newtons for force,
- pascals for Young's modulus and stress,
- kilograms per cubic meter for density.

External identifiers may be arbitrary `usize` values. Future solver stages will
map them to contiguous zero-based internal indices before assembly.

## Development Checks

The intended checks for every completed increment are:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```
