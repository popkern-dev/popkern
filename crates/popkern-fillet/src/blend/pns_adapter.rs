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
//!   `NurbsCurve` to a bi-2 Bézier — samples at `u = 0, 0.5, 1`
//!   and solves for the middle BB-coef so the Bézier interpolates
//!   those three points.
//! - **Derivative row** (3 CPs): the first transversal derivative
//!   of the supporting surface at the boundary — points INTO the
//!   surrounding surface (away from the cap), tangent to the
//!   surface, perpendicular to the boundary curve. Computed from
//!   the surface's partials `∂_u`, `∂_v` at the CP's closest-point
//!   UV, combined with the curve tangent.
//!
//! # Two entry points
//!
//! - [`blend_boundary_to_pns_ring`][]: placeholder derivatives
//!   (inward toward vertex). Fast path — no surface queries —
//!   useful for smoke tests. Cap is only G⁰ at the boundary.
//! - [`blend_boundary_to_pns_ring_with_surfaces`][]: surface-aware
//!   derivatives via `FilletSurfaceResult` + `BRep`. Produces a
//!   cap that matches supporting surface tangent planes at the
//!   boundary (G¹ target).
//!
//! # Valence
//!
//! If the blend has `N_fillets` fillets, `BlendBoundary` has
//! `2 · N_fillets` segments → the PnS cap has valence
//! `2·N_fillets`. Supported valences: `2·N_fillets ∈ {6, 8, 10, 12}`
//! (i.e., `N_fillets ∈ {3, 4, 5, 6}`).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::{Bi2BorderStrip, Bi2TensorBorderRing};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;

use super::boundary::{BlendBoundary, BoundarySegment};
use crate::fillet_surface::FilletSurfaceResult;

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
    /// A segment referenced a fillet or face whose supporting
    /// surface couldn't be accessed (out-of-range index, missing
    /// surface, etc.).
    #[error("PnS adapter: can't access supporting surface for segment {segment_idx}: {reason}")]
    SurfaceAccessFailed {
        /// Index of the failing segment.
        segment_idx: usize,
        /// Human-readable reason.
        reason: String,
    },
}

/// Fraction of segment length to use as the derivative magnitude.
/// Empirical default that gives a non-degenerate cap.
const DERIV_MAG_FRAC: f64 = 0.5;

/// Placeholder-derivative adapter. Faster (no surface queries) but
/// produces only a G⁰-matching cap. See module docs.
pub fn blend_boundary_to_pns_ring(
    boundary: &BlendBoundary,
) -> Result<Bi2TensorBorderRing, PnsAdapterError> {
    let n_seg = validate_segment_count(boundary.segments.len())?;
    let vertex = boundary.vertex_position;
    let mut sectors = Vec::with_capacity(n_seg);
    for segment in &boundary.segments {
        sectors.push(segment_to_strip_placeholder(segment, vertex));
    }
    Ok(Bi2TensorBorderRing { sectors })
}

/// Surface-aware adapter. Produces derivatives that encode the
/// supporting surface's tangent plane at each CP location → G¹
/// target when composed with the PnS cap assembly.
///
/// `fillet_surfaces` must be indexable by each FilletEnd segment's
/// `fillet_index`; `brep` must contain each FaceEdge segment's
/// `face` handle with a defined surface geometry.
pub fn blend_boundary_to_pns_ring_with_surfaces(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
) -> Result<Bi2TensorBorderRing, PnsAdapterError> {
    let n_seg = validate_segment_count(boundary.segments.len())?;
    let vertex = boundary.vertex_position;
    let mut sectors = Vec::with_capacity(n_seg);
    for (idx, segment) in boundary.segments.iter().enumerate() {
        sectors.push(segment_to_strip_surface_aware(
            segment,
            idx,
            fillet_surfaces,
            brep,
            vertex,
        )?);
    }
    Ok(Bi2TensorBorderRing { sectors })
}

fn validate_segment_count(n_seg: usize) -> Result<usize, PnsAdapterError> {
    if n_seg < 6 {
        return Err(PnsAdapterError::TooFewSegments(n_seg));
    }
    if n_seg > 12 {
        return Err(PnsAdapterError::TooManySegments(n_seg));
    }
    if !n_seg.is_multiple_of(2) {
        return Err(PnsAdapterError::OddSegmentCount(n_seg));
    }
    Ok(n_seg)
}

/// Sample each segment's curve at `t ∈ {t_min, t_mid, t_max}` and
/// return `(pos_bb, samples, curve_len)` where `pos_bb` is the bi-2
/// BB-coef row interpolating the 3 samples and `samples` is the
/// raw 3 sample points for downstream deriv computation.
fn pos_row_from_curve(segment: &BoundarySegment) -> ([Point3; 3], [Point3; 3], [f64; 3], f64) {
    let curve = match segment {
        BoundarySegment::FilletEnd { curve, .. } => curve,
        BoundarySegment::FaceEdge { curve, .. } => curve,
    };
    let (t_min, t_max) = curve.parameter_range();
    let t_mid = 0.5 * (t_min + t_max);
    let p_start = curve.point_at(t_min);
    let p_mid_sampled = curve.point_at(t_mid);
    let p_end = curve.point_at(t_max);
    // Bi-2 interpolating coef:  p_1 = 2·p(0.5) - 0.5·(p_0 + p_2)
    let p_mid_bb = Point3::new(
        2.0 * p_mid_sampled.x - 0.5 * (p_start.x + p_end.x),
        2.0 * p_mid_sampled.y - 0.5 * (p_start.y + p_end.y),
        2.0 * p_mid_sampled.z - 0.5 * (p_start.z + p_end.z),
    );
    let pos_bb = [p_start, p_mid_bb, p_end];
    let samples = [p_start, p_mid_sampled, p_end];
    let ts = [t_min, t_mid, t_max];
    let len = (p_end.x - p_start.x)
        .hypot(p_end.y - p_start.y)
        .hypot((p_end.z - p_start.z).abs());
    (pos_bb, samples, ts, len)
}

fn segment_to_strip_placeholder(segment: &BoundarySegment, vertex: Point3) -> Bi2BorderStrip {
    let (pos_bb, _samples, _ts, seg_len) = pos_row_from_curve(segment);
    let mag = seg_len * DERIV_MAG_FRAC;
    let deriv = [
        inward_deriv(pos_bb[0], vertex, mag),
        inward_deriv(pos_bb[1], vertex, mag),
        inward_deriv(pos_bb[2], vertex, mag),
    ];
    Bi2BorderStrip { pos: pos_bb, deriv }
}

fn segment_to_strip_surface_aware(
    segment: &BoundarySegment,
    segment_idx: usize,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    vertex: Point3,
) -> Result<Bi2BorderStrip, PnsAdapterError> {
    let surface = lookup_surface(segment, segment_idx, fillet_surfaces, brep)?;
    let curve = match segment {
        BoundarySegment::FilletEnd { curve, .. } => curve,
        BoundarySegment::FaceEdge { curve, .. } => curve,
    };
    let (pos_bb, samples, ts, seg_len) = pos_row_from_curve(segment);
    let mag = seg_len * DERIV_MAG_FRAC;

    let mut deriv = [Point3::new(0.0, 0.0, 0.0); 3];
    for i in 0..3 {
        let sample = samples[i];
        let Some((u, v)) = surface.closest_point(&sample) else {
            // Fallback: placeholder inward-pointing deriv.
            deriv[i] = inward_deriv(pos_bb[i], vertex, mag);
            continue;
        };
        let su = surface.partial_u(u, v);
        let sv = surface.partial_v(u, v);
        let n = su.cross(&sv);
        if n.length() < 1e-12 {
            deriv[i] = inward_deriv(pos_bb[i], vertex, mag);
            continue;
        }
        let curve_tangent = curve_tangent_at(curve, ts[i]);
        // In-surface perpendicular to curve tangent.
        let in_surf_perp = n.cross(&curve_tangent);
        let perp_len = in_surf_perp.length();
        if perp_len < 1e-12 {
            deriv[i] = inward_deriv(pos_bb[i], vertex, mag);
            continue;
        }
        let unit_perp = in_surf_perp * (1.0 / perp_len);
        // Orient to point AWAY from the blend vertex (into the
        // surrounding surface). If `unit_perp` currently points
        // toward the vertex, flip.
        let to_vertex = Vec3::new(
            vertex.x - sample.x,
            vertex.y - sample.y,
            vertex.z - sample.z,
        );
        let sign = if unit_perp.dot(&to_vertex) > 0.0 {
            -1.0
        } else {
            1.0
        };
        deriv[i] = Point3::new(
            sign * unit_perp.x * mag,
            sign * unit_perp.y * mag,
            sign * unit_perp.z * mag,
        );
    }
    Ok(Bi2BorderStrip { pos: pos_bb, deriv })
}

fn lookup_surface<'a>(
    segment: &'a BoundarySegment,
    segment_idx: usize,
    fillet_surfaces: &'a [FilletSurfaceResult],
    brep: &'a BRep,
) -> Result<&'a SurfaceGeometry, PnsAdapterError> {
    match segment {
        BoundarySegment::FilletEnd { fillet_index, .. } => fillet_surfaces
            .get(*fillet_index)
            .map(|fs| &fs.surface)
            .ok_or_else(|| PnsAdapterError::SurfaceAccessFailed {
                segment_idx,
                reason: format!("fillet_index {fillet_index} out of range"),
            }),
        BoundarySegment::FaceEdge { face, .. } => brep
            .faces
            .get(*face)
            .and_then(|f| f.surface.as_ref())
            .map(|sof| &sof.geometry)
            .ok_or_else(|| PnsAdapterError::SurfaceAccessFailed {
                segment_idx,
                reason: format!("face {face:?} has no surface geometry"),
            }),
    }
}

fn curve_tangent_at(curve: &popkern_core::geometry::NurbsCurve, t: f64) -> Vec3 {
    // First-order finite-difference tangent. The NurbsCurve trait
    // doesn't expose derivatives directly via CurveEval; sampling
    // two nearby points gives a robust enough approximation for
    // the adapter's purposes.
    let (t_min, t_max) = curve.parameter_range();
    let h = ((t_max - t_min) * 1e-5).max(1e-9);
    let t_lo = (t - h).max(t_min);
    let t_hi = (t + h).min(t_max);
    let p_lo = curve.point_at(t_lo);
    let p_hi = curve.point_at(t_hi);
    let dt = t_hi - t_lo;
    if dt.abs() < 1e-15 {
        return Vec3::new(1.0, 0.0, 0.0);
    }
    Vec3::new(
        (p_hi.x - p_lo.x) / dt,
        (p_hi.y - p_lo.y) / dt,
        (p_hi.z - p_lo.z) / dt,
    )
}

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
