// Stages 2b.31–2b.34 build infrastructure that's consumed when 2b.35
// wires the dispatcher into `try_dispatch_trig`. Until then, the
// inner functions look unused to the lib build. Remove this allow
// once 2b.35 lands.
#![allow(dead_code)]

//! Levin's algorithm for quadric × quadric intersection.
//!
//! Implements the closed-form parametric algorithm of Levin (1976/79)
//! "A parametric algorithm for drawing pictures of solid objects
//! composed of quadric surfaces", with the Wang-Goldman-Tu (2003)
//! robust treatment of degenerate pencils.
//!
//! ## Scope
//!
//! Closes the remaining analytic gaps in the EdgeTriple dispatch:
//! - Cylinder × cylinder skew (oblique non-parallel, non-perp axes)
//! - Cylinder × cone non-coplanar axes
//! - Cone × cone non-coaxial
//!
//! Specialised paths for parallel cylinders, axial cones, etc.
//! continue to take priority (they're exact and faster). Levin is the
//! fallback before NURBS marching.
//!
//! ## Pipeline
//!
//! 1. `quadric::surface_to_quadric` — `SurfaceGeometry` → `Mat4Sym`
//!    (the 4×4 symmetric matrix `A` such that `pᵀ·A·p = 0` on the
//!    surface, with `p = (x, y, z, w)` homogeneous).
//! 2. `pencil` — characteristic cubic `det(λA + μB)`, find a
//!    degenerate (rank ≤ 3) ruled quadric in the pencil.  *(stage
//!    2b.32)*
//! 3. `reduce` — parameterise the ruled quadric, substitute into one
//!    of the originals → quartic in one parameter.  *(stage 2b.33)*
//! 4. `solve` — Ferrari-Cardano on the quartic, classify each branch
//!    (ellipse / parabola / hyperbola) via plane-section
//!    discriminant.  *(stage 2b.34)*
//! 5. Dispatch wiring + cardinality on provenance.  *(stages 2b.35–37)*
//!
//! ## Symbolic-rational mode
//!
//! Stage 2b.36 lifts the f64 path to `Q[√d]` (rationals adjoined with
//! a single algebraic square root), giving exact branch classification
//! and parameterisation. Surd-bearing endpoints emit as f64 only at
//! the final pcurve / curve-3D layer.

mod quadric;

#[cfg(test)]
mod quadric_tests;
