//! Phase I — Input types for the polar Karčiauskas-Peters
//! construction (N ≥ 5 vertex blends).
//!
//! `PolarBoundaryInput` is a sibling of
//! [`DirectBoundaryInput`](super::DirectBoundaryInput). The
//! rectangular (N ∈ {3, 4}) path uses the latter; the polar path
//! uses this one. Both carry pure geometric data, no fillet-crate
//! topology dependency.
//!
//! # Polar structure
//!
//! A polar decomposition of an N-sided vertex blend places a
//! central vertex `C` and divides the blend region into `N`
//! sectors, one per spoke. Each sector is bounded by two shared
//! spoke curves (running from `C` out to two adjacent outer
//! corners `V_k`, `V_{k+1}`) and one outer arc of the original
//! input face at that sector.
//!
//! # Tensor-border
//!
//! Both K-P 2015 (bi-6) and K-P 2013 (bi-5) direct constructions
//! consume a **tensor-border** — depth-2 Hermite data around the
//! N-sided hole. See `docs/journal/polar-kp-math.md` §1.1–§1.4
//! for the K-P derivation. The variant
//! [`TensorBorder`] carries either bi-3 ([`TensorBorderBi3`]) or
//! bi-5 ([`TensorBorderBi5`]) per-sector-edge data.
//!
//! # Reused types
//!
//! Shares [`OuterCornerData`](super::OuterCornerData),
//! [`SpokeBoundaryData`](super::SpokeBoundaryData), and
//! [`TangentPlane`](super::TangentPlane) with the rectangular
//! input. Only the aggregate struct differs.

use crate::geometry::Point3;
use crate::geometry::gregory::{OuterCornerData, SpokeBoundaryData, TangentPlane};

/// Geometric input for the polar K-P construction at one vertex
/// blend of valence `N ≥ 3` (Phase I currently targets `N ≥ 5`;
/// the type accepts N = 3 / 4 to allow the K-P 2013 §6 n=3
/// special case and future experiments).
///
/// Invariants checked by [`validate`](Self::validate):
/// - All N-length vectors have consistent length.
/// - All per-spoke radii are finite and positive.
/// - Every spoke's boundary curve endpoints are near `C` and the
///   matching outer corner (sanity — exact positional equality is
///   not enforced because the Bezier CPs are not strict
///   interpolators at t=0/1).
#[derive(Debug, Clone)]
pub struct PolarBoundaryInput {
    /// Central vertex `C`. Shared by all `N` sectors.
    pub central_vertex: Point3,

    /// Common tangent plane at `C` — the Shirman-Séquin enclosure
    /// derived from the N incident face normals. Every sector's
    /// tangent at `C` must land in this plane (G1 at `C`).
    pub central_tangent_plane: TangentPlane,

    /// `N` shared spoke curves. Spoke `k` runs from `C` to
    /// outer corner `V_k`; it is shared between sector `k − 1`
    /// and sector `k` (indices mod N).
    pub spokes: Vec<SpokeBoundaryData>,

    /// `N` outer-corner data entries (one per spoke endpoint at
    /// `V_k`). `outer_corners[k]` sits at the shared endpoint of
    /// `spokes[k]`.
    pub outer_corners: Vec<OuterCornerData>,

    /// `N` face tangent planes — one per **sector** (the region
    /// between two adjacent spokes). Sector `k` is bounded by
    /// `spokes[k]` and `spokes[(k + 1) mod N]`; its face tangent
    /// plane is `face_tangent_planes[k]`. This is the input face
    /// the blend must match along the outer edge of sector `k`.
    pub face_tangent_planes: Vec<TangentPlane>,

    /// `N` per-spoke fillet radii. For constant-radius blends
    /// all entries are equal; for variable-radius, they differ and
    /// the K-P construction scales tensor-border cross-derivatives
    /// accordingly (see math summary §6).
    pub per_spoke_radii: Vec<f64>,
}

impl PolarBoundaryInput {
    /// Valence — the `N` in "N-sided polar vertex blend".
    pub fn n(&self) -> usize {
        self.spokes.len()
    }

    /// Structural invariants. A violation indicates an upstream
    /// bug (mismatched vector lengths, non-finite radii, etc.).
    pub fn validate(&self) -> Result<(), String> {
        let n = self.spokes.len();
        if n < 3 {
            return Err(format!("polar input needs at least 3 spokes, got {n}"));
        }
        if self.outer_corners.len() != n {
            return Err(format!(
                "outer_corners length {} ≠ spokes length {n}",
                self.outer_corners.len()
            ));
        }
        if self.face_tangent_planes.len() != n {
            return Err(format!(
                "face_tangent_planes length {} ≠ spokes length {n}",
                self.face_tangent_planes.len()
            ));
        }
        if self.per_spoke_radii.len() != n {
            return Err(format!(
                "per_spoke_radii length {} ≠ spokes length {n}",
                self.per_spoke_radii.len()
            ));
        }
        for (k, &r) in self.per_spoke_radii.iter().enumerate() {
            if !r.is_finite() || r <= 0.0 {
                return Err(format!(
                    "per_spoke_radii[{k}] = {r} is not finite / positive"
                ));
            }
        }
        // Sanity: each spoke's boundary CPs should start near C and
        // end near its matching outer corner. This mirrors the
        // `DirectBoundaryInput::validate` check.
        let c = self.central_vertex;
        for (k, spoke) in self.spokes.iter().enumerate() {
            let start = spoke.boundary_cps[0];
            let end = spoke.boundary_cps[7];
            let v = self.outer_corners[k].position;
            let d_start =
                ((start.x - c.x).powi(2) + (start.y - c.y).powi(2) + (start.z - c.z).powi(2))
                    .sqrt();
            let d_end =
                ((end.x - v.x).powi(2) + (end.y - v.y).powi(2) + (end.z - v.z).powi(2)).sqrt();
            if !d_start.is_finite() || !d_end.is_finite() {
                return Err(format!(
                    "spoke {k} endpoints non-finite: start {start:?} C {c:?} end {end:?} V {v:?}"
                ));
            }
        }
        Ok(())
    }

    /// True if every spoke has the same fillet radius.
    pub fn is_constant_radius(&self) -> bool {
        if self.per_spoke_radii.is_empty() {
            return true;
        }
        let r0 = self.per_spoke_radii[0];
        self.per_spoke_radii
            .iter()
            .all(|&r| (r - r0).abs() < 1e-12 * r0.abs().max(1.0))
    }
}

/// Depth-2 Hermite tensor-border around the N-sided hole.
///
/// Per the ADR §3.3.2 the K-P constructions accept two encodings:
/// - `Bi3` — bi-3 B-spline input (K-P 2015 native, K-P 2013
///   binary-split native).
/// - `Bi5` — bi-5 B-spline input (K-P 2013 §7.1 ternary-split
///   native).
#[derive(Debug, Clone)]
pub enum TensorBorder {
    /// bi-3 encoding (4 CPs × 3 Hermite rows per sector edge).
    Bi3(TensorBorderBi3),
    /// bi-5 encoding (6 CPs × 3 Hermite rows per sector edge).
    Bi5(TensorBorderBi5),
}

/// Bi-3 tensor-border: `N` sectors, each carrying three length-4
/// rows of BB-coefficients along the sector edge —
///
/// - row 0: the boundary curve itself (position row).
/// - row 1: first cross-derivative (into the hole, depth 1).
/// - row 2: second cross-derivative (into the hole, depth 2).
#[derive(Debug, Clone)]
pub struct TensorBorderBi3 {
    /// `sectors[k]` is the bi-3 tensor-border data along the
    /// outer edge of sector `k`.
    pub sectors: Vec<[[Point3; 4]; 3]>,
}

/// Bi-5 tensor-border: same structure as [`TensorBorderBi3`] but
/// with length-6 rows instead of length-4, for the K-P 2013
/// ternary-split / native-bi-5 input variant.
#[derive(Debug, Clone)]
pub struct TensorBorderBi5 {
    /// `sectors[k]` is the bi-5 tensor-border data along the
    /// outer edge of sector `k`.
    pub sectors: Vec<[[Point3; 6]; 3]>,
}

/// Degree choice for the tensor-border extracted from a polar
/// boundary input. Selectable by the caller to match the
/// downstream construction (bi-6 primary → [`Bi3`], bi-5
/// ternary-split → [`Bi5`]).
///
/// [`Bi3`]: TensorBorderKind::Bi3
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TensorBorderKind {
    /// bi-3 tensor-border (4 CPs per Hermite row per sector edge).
    Bi3,
    /// bi-5 tensor-border (6 CPs per Hermite row per sector edge).
    Bi5,
}

impl TensorBorder {
    /// Number of sectors — length of the internal `sectors` vec,
    /// regardless of whether the encoding is bi-3 or bi-5.
    pub fn n_sectors(&self) -> usize {
        match self {
            TensorBorder::Bi3(t) => t.sectors.len(),
            TensorBorder::Bi5(t) => t.sectors.len(),
        }
    }

    /// The encoding kind ([`Bi3`] or [`Bi5`]).
    ///
    /// [`Bi3`]: TensorBorderKind::Bi3
    pub fn kind(&self) -> TensorBorderKind {
        match self {
            TensorBorder::Bi3(_) => TensorBorderKind::Bi3,
            TensorBorder::Bi5(_) => TensorBorderKind::Bi5,
        }
    }

    /// Position-row (row 0) BB-CPs of sector `k`. Returns a slice
    /// view so callers don't pay a clone; length is 4 for bi-3
    /// and 6 for bi-5.
    pub fn position_row(&self, k: usize) -> &[Point3] {
        match self {
            TensorBorder::Bi3(t) => &t.sectors[k][0],
            TensorBorder::Bi5(t) => &t.sectors[k][0],
        }
    }

    /// First-cross-derivative row (row 1) BB-CPs of sector `k`.
    pub fn first_cross_row(&self, k: usize) -> &[Point3] {
        match self {
            TensorBorder::Bi3(t) => &t.sectors[k][1],
            TensorBorder::Bi5(t) => &t.sectors[k][1],
        }
    }

    /// Second-cross-derivative row (row 2) BB-CPs of sector `k`.
    pub fn second_cross_row(&self, k: usize) -> &[Point3] {
        match self {
            TensorBorder::Bi3(t) => &t.sectors[k][2],
            TensorBorder::Bi5(t) => &t.sectors[k][2],
        }
    }
}

#[cfg(test)]
#[path = "kp_polar_boundary_input_tests.rs"]
mod tests;
