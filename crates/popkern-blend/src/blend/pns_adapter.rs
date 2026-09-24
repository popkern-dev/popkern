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
//!   `∂_v b(u, 0)` of the tensor border, where `v` increases INTO
//!   the cap: tangent to the supporting surface, perpendicular to
//!   the boundary curve, pointing across the boundary toward the
//!   cap interior (the direction the surrounding surface continues
//!   into the hole). Computed from the surface's partials `∂_u`,
//!   `∂_v` at the CP's closest-point UV, combined with the curve
//!   tangent. The cap consumes this row two ways, both of which
//!   need the inward sense: `b ∘ β` composes the cap's first
//!   interior row as `b·∂_u b + a·∂_v b` with `a(u) > 0`, and the
//!   pole (eq 5) is a convex combination of `pos + deriv/2` — the
//!   net's inner ring. An outward row (the pre-PK-517 convention)
//!   folds the cap's first row back over the fillet and drags the
//!   pole below the rim on every convex corner (measured on the
//!   carved-pentagon boundary ring: pole z 0.9528 under a rim at
//!   0.9603; inward: 0.9995, above the rim's arc crowns at 0.9762).
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
//! The adapter pairs each `(FilletEnd arc, FaceEdge bridge)` pair
//! into ONE PnS sector, so the PnS cap has valence `N_fillets`
//! (not `2·N_fillets`). This matches the `N`-sub-face output of
//! `poke_n_gon_face_into_n_quads` on a `2N`-HE blend face.
//!
//! Accepted valences: `N_fillets ∈ [3, 12]`, matching the cap
//! assembly (`n = 4`'s σ pole is removable — see
//! `popkern_core::geometry::gregory::pns_sigma`). Acceptance here is
//! *representational*; which valences the production dispatch
//! actually serves is decided by measurement. PK-192's recorded
//! route comparison keeps `n ∈ {3, 4}` on the direct-apex/Coons
//! route: pairing arc + bridge into one sector means each border
//! row must fit a kinked curve pair with a single degree-4 row, a
//! structural floor of ≈ 6e-3 against the 5e-3 acceptance bar
//! (pyramid3 6.458e-3, cube corner 5.972e-3, support exact), while
//! the incumbent interpolates the same rim at ~1e-15. See
//! `pns_low_valence_gates.rs` for the decision record and the
//! reversal condition.

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
    /// Too many segments — PnS valence capped at 12 (24 segments).
    #[error("PnS adapter: boundary has {0} segments, expected ≤ 24 (valence ≤ 12)")]
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
///
/// Pairs consecutive segments `(2k, 2k+1)` into one PnS sector.
pub fn blend_boundary_to_pns_ring(
    boundary: &BlendBoundary,
) -> Result<Bi2TensorBorderRing, PnsAdapterError> {
    let n = validate_and_compute_valence(boundary.segments.len())?;
    let vertex = boundary.vertex_position;
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let seg_a = &boundary.segments[2 * k];
        let seg_b = &boundary.segments[2 * k + 1];
        sectors.push(pair_to_strip_placeholder(seg_a, seg_b, vertex));
    }
    Ok(Bi2TensorBorderRing { sectors })
}

/// Surface-aware adapter. Produces derivatives that encode the
/// supporting surface's tangent plane → G¹ target.
///
/// Pairs consecutive segments `(2k, 2k+1)` into one PnS sector:
/// `pos[0]` at arc start, `pos[1]` at arc/bridge join, `pos[2]` at
/// bridge end. Derivatives `deriv[0]`, `deriv[1]` from the arc's
/// surface; `deriv[2]` from the bridge's surface.
pub fn blend_boundary_to_pns_ring_with_surfaces(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
) -> Result<Bi2TensorBorderRing, PnsAdapterError> {
    let n = validate_and_compute_valence(boundary.segments.len())?;
    let vertex = boundary.vertex_position;
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let seg_a = &boundary.segments[2 * k];
        let seg_b = &boundary.segments[2 * k + 1];
        sectors.push(pair_to_strip_surface_aware(
            seg_a,
            seg_b,
            2 * k,
            fillet_surfaces,
            brep,
            vertex,
        )?);
    }
    Ok(Bi2TensorBorderRing { sectors })
}

/// Validate segment count and return the PnS valence
/// `n = n_seg / 2`.
fn validate_and_compute_valence(n_seg: usize) -> Result<usize, PnsAdapterError> {
    if n_seg < 6 {
        return Err(PnsAdapterError::TooFewSegments(n_seg));
    }
    if n_seg > 24 {
        return Err(PnsAdapterError::TooManySegments(n_seg));
    }
    if !n_seg.is_multiple_of(2) {
        return Err(PnsAdapterError::OddSegmentCount(n_seg));
    }
    Ok(n_seg / 2)
}

/// Build a bi-2 position row `[p_start_a, p_join, p_end_b]` that
/// interpolates the two segments' endpoints AND the arc/bridge
/// join point at parameter `u = 0.5`.
fn paired_pos_row(seg_a: &BoundarySegment, seg_b: &BoundarySegment) -> [Point3; 3] {
    let curve_a = curve_of(seg_a);
    let curve_b = curve_of(seg_b);
    let (a_min, a_max) = curve_a.parameter_range();
    let (_, b_max) = curve_b.parameter_range();
    let p_start = curve_a.point_at(a_min);
    let p_join = curve_a.point_at(a_max); // = curve_b.point_at(b_min) by cyclic consistency
    let p_end = curve_b.point_at(b_max);
    // Bi-2 interpolating coef at u=0.5:
    //   p_join = 0.25·p_start + 0.5·p_mid_bb + 0.25·p_end
    // → p_mid_bb = 2·p_join − 0.5·(p_start + p_end)
    let p_mid_bb = Point3::new(
        2.0 * p_join.x - 0.5 * (p_start.x + p_end.x),
        2.0 * p_join.y - 0.5 * (p_start.y + p_end.y),
        2.0 * p_join.z - 0.5 * (p_start.z + p_end.z),
    );
    [p_start, p_mid_bb, p_end]
}

fn curve_of(segment: &BoundarySegment) -> &popkern_core::geometry::NurbsCurve {
    match segment {
        BoundarySegment::FilletEnd { curve, .. } => curve,
        BoundarySegment::FaceEdge { curve, .. } => curve,
    }
}

fn pair_to_strip_placeholder(
    seg_a: &BoundarySegment,
    seg_b: &BoundarySegment,
    vertex: Point3,
) -> Bi2BorderStrip {
    let pos_bb = paired_pos_row(seg_a, seg_b);
    // Rough magnitude: combined arc+bridge length × fraction.
    let curve_a = curve_of(seg_a);
    let curve_b = curve_of(seg_b);
    let (a_min, a_max) = curve_a.parameter_range();
    let (b_min, b_max) = curve_b.parameter_range();
    let len_a = curve_a
        .point_at(a_min)
        .distance_to(&curve_a.point_at(a_max));
    let len_b = curve_b
        .point_at(b_min)
        .distance_to(&curve_b.point_at(b_max));
    let mag = (len_a + len_b) * DERIV_MAG_FRAC;
    let deriv = [
        inward_deriv(pos_bb[0], vertex, mag),
        inward_deriv(pos_bb[1], vertex, mag),
        inward_deriv(pos_bb[2], vertex, mag),
    ];
    Bi2BorderStrip { pos: pos_bb, deriv }
}

fn pair_to_strip_surface_aware(
    seg_a: &BoundarySegment,
    seg_b: &BoundarySegment,
    base_idx: usize,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    vertex: Point3,
) -> Result<Bi2BorderStrip, PnsAdapterError> {
    // Soft lookup: if either segment's supporting surface isn't
    // accessible (common in mid-pipeline BReps where internal
    // faces don't yet have surface geometry attached), fall back
    // to placeholder derivatives for that segment's samples.
    let surface_a = lookup_surface(seg_a, base_idx, fillet_surfaces, brep).ok();
    let surface_b = lookup_surface(seg_b, base_idx + 1, fillet_surfaces, brep).ok();
    let curve_a = curve_of(seg_a);
    let curve_b = curve_of(seg_b);
    let (a_min, a_max) = curve_a.parameter_range();
    let (b_min, b_max) = curve_b.parameter_range();

    let pos_bb = paired_pos_row(seg_a, seg_b);
    let len_a = curve_a
        .point_at(a_min)
        .distance_to(&curve_a.point_at(a_max));
    let len_b = curve_b
        .point_at(b_min)
        .distance_to(&curve_b.point_at(b_max));
    let mag = (len_a + len_b) * DERIV_MAG_FRAC;

    // Per-position optional surface + sample-t + associated curve.
    let sample_data: [(
        Option<&popkern_core::geometry::surface::SurfaceGeometry>,
        &popkern_core::geometry::NurbsCurve,
        f64,
        Point3,
    ); 3] = [
        (surface_a, curve_a, a_min, pos_bb[0]),
        (surface_a, curve_a, a_max, pos_bb[1]),
        (surface_b, curve_b, b_max, pos_bb[2]),
    ];

    let mut deriv = [Point3::new(0.0, 0.0, 0.0); 3];
    for (i, (surface_opt, curve, t, sample_pt)) in sample_data.iter().enumerate() {
        let Some(surface) = surface_opt else {
            deriv[i] = inward_deriv(pos_bb[i], vertex, mag);
            continue;
        };
        let Some((u, v)) = surface.closest_point(sample_pt) else {
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
        let curve_tangent = curve_tangent_at(curve, *t);
        let in_surf_perp = n.cross(&curve_tangent);
        let perp_len = in_surf_perp.length();
        if perp_len < 1e-12 {
            deriv[i] = inward_deriv(pos_bb[i], vertex, mag);
            continue;
        }
        let unit_perp = in_surf_perp * (1.0 / perp_len);
        let to_vertex = Vec3::new(
            vertex.x - sample_pt.x,
            vertex.y - sample_pt.y,
            vertex.z - sample_pt.z,
        );
        // Inward (toward the cap): the cap interior lies on the
        // vertex side of every boundary segment of a corner blend.
        let sign = if unit_perp.dot(&to_vertex) > 0.0 {
            1.0
        } else {
            -1.0
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
