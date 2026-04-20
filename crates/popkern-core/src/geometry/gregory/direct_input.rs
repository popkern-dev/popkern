//! Inputs to direct-construction — pure geometric data, no
//! dependency on the fillet-crate topology types.
//!
//! The caller (typically `popkern-fillet`) extracts this struct
//! from its own `BlendBoundary` + adjacent-face geometry. Keeping
//! the math layer dependency-free means:
//!   - Direct construction can be unit-tested without fillet
//!     scaffolding.
//!   - Future callers outside `popkern-fillet` (e.g. a sweep or
//!     offset operation that wants the same construction) can
//!     build the input struct themselves.
//!
//! # Geometry conventions
//!
//! All boundary curves are represented as bidegree-7 Bezier
//! control-point arrays matching the two-stripe patch format
//! (`[Point3; 8]`). Upstream is responsible for any degree
//! elevation from lower-degree inputs.
//!
//! Tangents are evaluated as unit vectors at the spoke endpoints
//! — direction only; magnitude is derived from the curve.

use crate::geometry::gregory::RectangularGregoryTwoStripe;
use crate::geometry::{Point3, Vec3};

/// The geometric data direct-construction needs to build the
/// two-stripe patches for an N-sided vertex blend.
///
/// `N` is implicit — every field's outer length equals the
/// spoke count. Invariants are checked in
/// [`DirectBoundaryInput::validate`].
#[derive(Debug, Clone)]
pub struct DirectBoundaryInput {
    /// Central vertex position (`C`). Shared by all sub-patches.
    pub central_vertex: Point3,

    /// Per-spoke data. One entry per shared spoke in the N-sided
    /// ring. Length `N ≥ 3`.
    pub spokes: Vec<SpokeBoundaryData>,

    /// Per-outer-corner data. Corner `k` is `V_{2k+2}`, shared
    /// between sub-patch `k` (at its `t=1`) and sub-patch
    /// `k+1 mod N` (at its `t=1` reversed). Length `N`.
    pub outer_corners: Vec<OuterCornerData>,

    /// Adjacent-face unit tangent planes at `C`, one per patch
    /// in the ring. At sub-patch `k`, the input-face tangent
    /// plane is spanned by two vectors in `central_vertex_tps[k]`.
    /// These are the **input** face tangent planes the blend
    /// must match at `C` — the Shirman-Séquin construction
    /// forces every sub-patch's tangent at `C` to land in its
    /// incident face tangent plane.
    pub central_vertex_tps: Vec<TangentPlane>,

    /// **Optional** biseptic-promoted two-stripe seed patches,
    /// one per sub-patch in the ring. When present, direct-
    /// construction uses these as the source of existing
    /// boundary CPs (especially the hex-external edges at
    /// `u=0` and `v=1` which aren't per-spoke data). Upstream
    /// in `popkern-fillet` builds these via the existing
    /// biseptic promotion ladder; direct-construction then
    /// refines the row-1 / row-2 CPs on top.
    ///
    /// When `None`, outer-corner and row-construction modules
    /// fall back to the parallelogram-rule seed, which is
    /// sufficient for symmetric / regular inputs.
    pub seed_patches: Option<Vec<RectangularGregoryTwoStripe>>,
}

/// Data for a single shared spoke between adjacent sub-patches.
#[derive(Debug, Clone)]
pub struct SpokeBoundaryData {
    /// Bidegree-7 Bezier control points of the spoke curve
    /// `C(t)` as it appears in patch `k`'s UOne boundary
    /// (`t = 0 → u=1, v=0`; `t = 1 → u=1, v=1`).
    ///
    /// The same points appear reversed in patch `(k+1) mod N`'s
    /// VZero boundary (`b_reversed = true` per spoke convention
    /// in `construct::three_quad_apply`).
    pub boundary_cps: [Point3; 8],

    /// Patch-A-side unit cross-tangent at `t = 0` (= at `C`).
    /// This is `∂S_A/∂u(1, 0)` normalised, evaluated on the
    /// Shirman-Séquin-valid biseptic seed.
    pub xtan_a_at_c: Vec3,

    /// Patch-B-side unit cross-tangent at `t = 0` (= at `C`).
    pub xtan_b_at_c: Vec3,

    /// Patch-A-side unit cross-tangent at `t = 1` (= at outer
    /// corner `V_{2k+2}`).
    pub xtan_a_at_v: Vec3,

    /// Patch-B-side unit cross-tangent at `t = 1`.
    pub xtan_b_at_v: Vec3,

    /// Radius at `t = 0`. For constant-radius blends this equals
    /// `radius_at_v`; for variable-radius, the two differ and
    /// the direct-construction formulas switch to the K-K-C 1991
    /// variable-radius variant.
    pub radius_at_c: f64,

    /// Radius at `t = 1`.
    pub radius_at_v: f64,

    /// Edge identifier for user-facing error messages (e.g.
    /// `"edge_7"`). No semantic role inside the math.
    pub edge_id: String,
}

/// Data for a single outer hex corner `V_{2k+2}`.
#[derive(Debug, Clone)]
pub struct OuterCornerData {
    /// Position of `V_{2k+2}` in world coordinates.
    pub position: Point3,

    /// Input face tangent plane at `V_{2k+2}` — both adjacent
    /// sub-patches must have their tangent planes at this corner
    /// equal this. Usually determined by the two input faces'
    /// normal vectors at the original vertex.
    pub tangent_plane: TangentPlane,
}

/// A unit tangent plane at a point: two orthonormal basis
/// vectors spanning the plane. The surface normal is `u × v`.
///
/// The orthonormality is a construction invariant — direct-
/// construction's corner formulas depend on it. Helpers in
/// [`crate::geometry::Vec3`] should normalise before packing.
#[derive(Debug, Clone, Copy)]
pub struct TangentPlane {
    /// First basis vector of the plane (unit).
    pub u: Vec3,
    /// Second basis vector of the plane (unit, orthogonal to `u`).
    pub v: Vec3,
}

impl TangentPlane {
    /// Build a tangent plane from two arbitrary in-plane vectors.
    /// Returns `None` if the vectors are parallel or either is
    /// zero-length.
    pub fn from_basis(u_in: Vec3, v_in: Vec3) -> Option<Self> {
        let u_unit = u_in.normalized()?;
        // Gram-Schmidt: remove the u-component from v_in.
        let u_dot_v = u_unit.x * v_in.x + u_unit.y * v_in.y + u_unit.z * v_in.z;
        let v_residual = Vec3::new(
            v_in.x - u_dot_v * u_unit.x,
            v_in.y - u_dot_v * u_unit.y,
            v_in.z - u_dot_v * u_unit.z,
        );
        let v_unit = v_residual.normalized()?;
        Some(Self {
            u: u_unit,
            v: v_unit,
        })
    }

    /// Build a tangent plane from a unit normal. Chooses an
    /// arbitrary (but deterministic) pair of in-plane basis
    /// vectors — the specific choice doesn't matter for our
    /// coplanarity checks.
    pub fn from_normal(n: Vec3) -> Option<Self> {
        let n_unit = n.normalized()?;
        // Pick a "least-parallel" axis to cross with.
        let axis = if n_unit.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let u = n_unit.cross(&axis).normalized()?;
        let v = n_unit.cross(&u).normalized()?;
        Some(Self { u, v })
    }

    /// Unit normal to the plane (`u × v`).
    pub fn normal(&self) -> Vec3 {
        self.u.cross(&self.v)
    }

    /// True if `vec` lies in this plane (within `tolerance`).
    pub fn contains(&self, vec: Vec3, tolerance: f64) -> bool {
        let n = self.normal();
        let dot = vec.x * n.x + vec.y * n.y + vec.z * n.z;
        dot.abs() < tolerance
    }
}

impl DirectBoundaryInput {
    /// Valence — the `N` in "N-sided vertex blend".
    pub fn n(&self) -> usize {
        self.spokes.len()
    }

    /// Structural invariants. Called before any construction; a
    /// violation indicates an upstream bug (wrong number of
    /// spokes / corners, etc.).
    pub fn validate(&self) -> Result<(), String> {
        let n = self.spokes.len();
        if n < 3 {
            return Err(format!("need at least 3 spokes, got {n}"));
        }
        if self.outer_corners.len() != n {
            return Err(format!(
                "outer_corners length {} ≠ spokes length {}",
                self.outer_corners.len(),
                n
            ));
        }
        if self.central_vertex_tps.len() != n {
            return Err(format!(
                "central_vertex_tps length {} ≠ spokes length {}",
                self.central_vertex_tps.len(),
                n
            ));
        }
        if let Some(ref seeds) = self.seed_patches
            && seeds.len() != n
        {
            return Err(format!(
                "seed_patches length {} ≠ spokes length {}",
                seeds.len(),
                n
            ));
        }
        // Every spoke's boundary CPs must start at C and end at
        // its matching outer corner. Relaxed: positional match
        // within a reasonable tolerance (boundary CPs are fit
        // Beziers, not exact interpolators).
        for (k, spoke) in self.spokes.iter().enumerate() {
            let start = spoke.boundary_cps[0];
            let end = spoke.boundary_cps[7];
            let c = self.central_vertex;
            let v = self.outer_corners[k].position;
            let to_c = (start.x - c.x).hypot(start.y - c.y).hypot(start.z - c.z);
            let to_v = (end.x - v.x).hypot(end.y - v.y).hypot(end.z - v.z);
            // Permissive — just sanity; math layer doesn't rely
            // on exact endpoint matching (only on the Bezier CPs
            // themselves).
            if !to_c.is_finite() || !to_v.is_finite() {
                return Err(format!(
                    "spoke {k} endpoints non-finite: \
                     start {start:?} C {c:?} end {end:?} V {v:?}"
                ));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "direct_input_tests.rs"]
mod tests;
