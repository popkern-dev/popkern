//! Adapter: `BlendBoundary` → PnS `Bi2TensorBorderRing`.
//!
//! Converts the CAD-level blend boundary (2N alternating
//! FilletEnd / FaceEdge segments around a vertex) into the PnS
//! cap's bi-2 tensor-border input (2N sectors).
//!
//! # Mapping
//!
//! Each `BoundarySegment` becomes one `Bi2BorderStrip`:
//! - **Position row** (3 CPs): degree-reduction of the segment's
//!   `NurbsCurve` to a bi-2 (degree-2) Bézier — samples at the
//!   curve's parameter endpoints and midpoint, then solves for
//!   the middle BB-coef so the Bézier interpolates those three
//!   points.
//! - **Derivative row** (3 CPs): **placeholder** — an inward-
//!   pointing vector from each position CP toward the blend
//!   vertex, scaled by a constant fraction of the segment length.
//!   This gives a sensible non-zero transverse direction without
//!   requiring access to the supporting surface (which would be
//!   needed for proper G¹ matching).
//!
//! **G¹ caveat**: with placeholder derivatives, the resulting
//! cap matches the blend boundary positions but NOT the surface
//! tangent planes. Proper derivatives require plumbing each
//! segment's supporting surface normal — deferred until we have a
//! concrete use-case for end-to-end G¹.
//!
//! # Valence
//!
//! If the blend has `N_fillets` fillets, `BlendBoundary` has
//! `2 · N_fillets` segments → the PnS cap has valence `2·N_fillets`.
//! Supported valences: `2·N_fillets ∈ {6, 8, 10, 12}` (i.e.,
//! `N_fillets ∈ {3, 4, 5, 6}`). `2·N = 4` (a 2-fillet blend) is
//! rejected by the downstream PnS solver — saddle vertices aren't
//! the common case.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::{Bi2BorderStrip, Bi2TensorBorderRing};

use super::boundary::{BlendBoundary, BoundarySegment};

/// Error cases for the `BlendBoundary → Bi2TensorBorderRing`
/// adapter.
#[derive(Debug, thiserror::Error)]
pub enum PnsAdapterError {
    /// Too few segments — PnS needs valence ≥ 3, i.e., ≥ 6
    /// alternating segments.
    #[error("PnS adapter: boundary has {0} segments, expected ≥ 6 (valence ≥ 3)")]
    TooFewSegments(usize),
    /// Too many segments — PnS valence capped at 6 (12 segments).
    #[error("PnS adapter: boundary has {0} segments, expected ≤ 12 (valence ≤ 6)")]
    TooManySegments(usize),
    /// Segment count is odd — `BlendBoundary` should have an even
    /// count (alternating FilletEnd/FaceEdge).
    #[error(
        "PnS adapter: boundary has {0} segments, expected even (alternating FilletEnd/FaceEdge)"
    )]
    OddSegmentCount(usize),
}

/// Adapter scaling factor for the placeholder transverse derivative.
/// Determines how far "into the blend" the derivative points, as a
/// fraction of each segment's arc length. `0.5` = half a segment
/// length — a reasonable default that produces non-degenerate
/// transverse data.
const DERIV_INWARD_FRAC: f64 = 0.5;

/// Convert a `BlendBoundary` to a bi-2 tensor border ring suitable
/// for `assemble_pns_bi4_cap`.
///
/// The resulting ring has `N = boundary.segments.len()` sectors,
/// one per `BoundarySegment`. Segment `k` becomes sector `k`'s
/// `Bi2BorderStrip`.
///
/// See the module docs for the mapping semantics and G¹ caveat.
pub fn blend_boundary_to_pns_ring(
    boundary: &BlendBoundary,
) -> Result<Bi2TensorBorderRing, PnsAdapterError> {
    let n_seg = boundary.segments.len();
    if n_seg < 6 {
        return Err(PnsAdapterError::TooFewSegments(n_seg));
    }
    if n_seg > 12 {
        return Err(PnsAdapterError::TooManySegments(n_seg));
    }
    if !n_seg.is_multiple_of(2) {
        return Err(PnsAdapterError::OddSegmentCount(n_seg));
    }

    let vertex = boundary.vertex_position;
    let mut sectors = Vec::with_capacity(n_seg);
    for segment in &boundary.segments {
        sectors.push(segment_to_strip(segment, vertex));
    }
    Ok(Bi2TensorBorderRing { sectors })
}

/// Build one `Bi2BorderStrip` from one `BoundarySegment`.
fn segment_to_strip(segment: &BoundarySegment, vertex: Point3) -> Bi2BorderStrip {
    let curve = match segment {
        BoundarySegment::FilletEnd { curve, .. } => curve,
        BoundarySegment::FaceEdge { curve, .. } => curve,
    };
    let (t_min, t_max) = curve.parameter_range();
    let t_mid = 0.5 * (t_min + t_max);

    let p_start = curve.point_at(t_min);
    let p_end = curve.point_at(t_max);
    let p_mid_sampled = curve.point_at(t_mid);

    // Solve for the bi-2 middle BB-coef so the Bézier
    //   p(u) = (1-u)² p_0 + 2u(1-u) p_1 + u² p_2
    // interpolates (p_start, p_mid_sampled, p_end) at u = 0, 0.5, 1.
    // At u=0.5: p = 0.25·p_0 + 0.5·p_1 + 0.25·p_2
    //        → p_1 = 2·p(0.5) - 0.5·(p_0 + p_2)
    let p_mid_bb = Point3::new(
        2.0 * p_mid_sampled.x - 0.5 * (p_start.x + p_end.x),
        2.0 * p_mid_sampled.y - 0.5 * (p_start.y + p_end.y),
        2.0 * p_mid_sampled.z - 0.5 * (p_start.z + p_end.z),
    );

    let pos = [p_start, p_mid_bb, p_end];

    // Placeholder deriv: each position CP gets an inward-pointing
    // vector (toward vertex), magnitude = segment length × fraction.
    let seg_len = (p_end.x - p_start.x)
        .hypot(p_end.y - p_start.y)
        .hypot((p_end.z - p_start.z).abs());
    let deriv_mag = seg_len * DERIV_INWARD_FRAC;
    let deriv = [
        inward_deriv(pos[0], vertex, deriv_mag),
        inward_deriv(pos[1], vertex, deriv_mag),
        inward_deriv(pos[2], vertex, deriv_mag),
    ];

    Bi2BorderStrip { pos, deriv }
}

/// Unit vector from `p` toward `vertex`, scaled by `mag`. Returns
/// zero if `p ≈ vertex`.
fn inward_deriv(p: Point3, vertex: Point3, mag: f64) -> Point3 {
    let dx = vertex.x - p.x;
    let dy = vertex.y - p.y;
    let dz = vertex.z - p.z;
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    if len < 1e-12 {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let s = mag / len;
    Point3::new(dx * s, dy * s, dz * s)
}
