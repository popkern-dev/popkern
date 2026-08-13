//! Feasibility pre-flight helpers extracted from
//! `geometry_compute.rs` to keep the parent under the 500-line
//! structural cap.
//!
//! Contains the self-intersection curvature check that the fillet
//! engine runs before attempting any real geometry work — if the
//! requested radius would make the offset of either support surface
//! self-intersect, we short-circuit with a structured
//! `BlendFailure::RadiusTooLarge` carrying the tighter bound so the
//! caller can reduce the radius and retry.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle};

use crate::offset::self_intersect::check_offset_self_intersection;
use crate::types::FilletRadius;

/// [`self_intersection_safe_radius`] with both face domains restricted
/// to the edge's contact band — the form the pipeline pre-flight and
/// the §10.3 precheck both use.
#[allow(clippy::too_many_arguments)]
pub(crate) fn banded_self_intersection_safe_radius(
    brep: &BRep,
    edge: EdgeHandle,
    surface_1: &SurfaceGeometry,
    surface_2: &SurfaceGeometry,
    radius: &FilletRadius,
    orientation_1: bool,
    orientation_2: bool,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> Option<f64> {
    let peak = largest_requested_radius(radius);
    let (u1, v1) = edge_contact_band(brep, edge, surface_1, u1_range, v1_range, peak);
    let (u2, v2) = edge_contact_band(brep, edge, surface_2, u2_range, v2_range, peak);
    self_intersection_safe_radius(
        surface_1,
        surface_2,
        radius,
        orientation_1,
        orientation_2,
        u1,
        v1,
        u2,
        v2,
    )
}

/// Restrict a support face's UV domain to the band the rolling ball
/// actually contacts around the filleted edge (PK-264).
///
/// The offset self-intersection pre-flight is a LOCAL curvature
/// feasibility question — "can a ball of this radius roll along this
/// edge?" — but it used to be evaluated over the face's entire
/// parameter domain, so any remote feature of the face capped the
/// radius for every edge on it. Measured: the cone base-circle fillet
/// declined RadiusTooLarge at max 0.0369 on a radius-1.0 rim because
/// the INWARD offset cone self-intersects near the APEX — a region a
/// base fillet never touches.
///
/// The band is the edge's own UV bounding box (branch-normalized
/// projection) padded by the chart-metric equivalent of `3 ×
/// peak_radius` — the contact strip sits within ~1 radius of the
/// edge, so 3× is generous — then intersected with the face domain.
/// Any failure to sample or project falls back to the full face
/// domain (the conservative pre-fix behavior).
pub(crate) fn edge_contact_band(
    brep: &BRep,
    edge: EdgeHandle,
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    peak_radius: f64,
) -> ((f64, f64), (f64, f64)) {
    use popkern_core::geometry::curve::CurveEval;
    let full = (u_range, v_range);
    let Some(curve) = brep.edges.get(edge).and_then(|e| e.curve.as_ref()) else {
        return full;
    };
    const SAMPLES: usize = 9;
    let mut u_lo = f64::INFINITY;
    let mut u_hi = f64::NEG_INFINITY;
    let mut v_lo = f64::INFINITY;
    let mut v_hi = f64::NEG_INFINITY;
    let mut hits = 0usize;
    for i in 0..SAMPLES {
        let t = curve.t_start + (curve.t_end - curve.t_start) * (i as f64) / ((SAMPLES - 1) as f64);
        let p = curve.geometry.point_at(t);
        if let Some((u, v)) =
            crate::center_curve_helpers::project_to_face_branch(surface, u_range, v_range, &p)
        {
            u_lo = u_lo.min(u);
            u_hi = u_hi.max(u);
            v_lo = v_lo.min(v);
            v_hi = v_hi.max(v);
            hits += 1;
        }
    }
    if hits < 2 {
        return full;
    }
    // Chart-metric pad: convert 3·r of arc length into parameter
    // units via the local derivative magnitudes at the band center.
    let (uc, vc) = (0.5 * (u_lo + u_hi), 0.5 * (v_lo + v_hi));
    let su = surface.partial_u(uc, vc).length();
    let sv = surface.partial_v(uc, vc).length();
    let pad = 3.0 * peak_radius;
    let pad_u = if su > 1.0e-12 { pad / su } else { pad };
    let pad_v = if sv > 1.0e-12 { pad / sv } else { pad };
    let band_u = ((u_lo - pad_u).max(u_range.0), (u_hi + pad_u).min(u_range.1));
    let band_v = ((v_lo - pad_v).max(v_range.0), (v_hi + pad_v).min(v_range.1));
    if band_u.0 >= band_u.1 || band_v.0 >= band_v.1 {
        return full;
    }
    (band_u, band_v)
}

/// Largest requested radius in a radius spec. Used as the "requested"
/// value in `BlendFailure::RadiusTooLarge`, so the user sees the
/// worst-case instance of their spec, not an arbitrary midpoint.
pub(crate) fn largest_requested_radius(radius: &FilletRadius) -> f64 {
    match radius {
        FilletRadius::Constant(r) => *r,
        FilletRadius::Variable { start, end } => start.max(*end),
        FilletRadius::ByPoints(pts) => pts
            .iter()
            .map(|(_, r)| *r)
            .fold(0.0_f64, |acc, r| acc.max(r)),
    }
}

/// Pre-flight self-intersection check for a rolling-ball fillet at the
/// requested radius. Returns `Some(max_safe_radius)` when the
/// requested radius would exceed the feasibility envelope on either
/// support surface, or `None` when the radius is safe (i.e. the
/// caller should proceed with the normal geometry pipeline).
///
/// For [`FilletRadius::Constant`] we test the single value. For
/// variable/by-points specs we test the peak requested radius —
/// if the peak is unsafe, the spec is unconditionally unsafe.
#[allow(clippy::too_many_arguments)]
pub(crate) fn self_intersection_safe_radius(
    surface_1: &SurfaceGeometry,
    surface_2: &SurfaceGeometry,
    radius: &FilletRadius,
    orientation_1: bool,
    orientation_2: bool,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> Option<f64> {
    let peak = largest_requested_radius(radius);
    if peak <= 0.0 {
        return None;
    }

    // "inward" matches the logic in center_curve.rs — offset_surface's
    // inward flag moves AGAINST the normal, and face normals point
    // into the solid, so we need `inward = !orientation`.
    let inward_1 = !orientation_1;
    let inward_2 = !orientation_2;

    let max_1 = check_offset_self_intersection(surface_1, peak, inward_1, u1_range, v1_range);
    let max_2 = check_offset_self_intersection(surface_2, peak, inward_2, u2_range, v2_range);

    // Both surfaces reported "safe" — nothing to flag.
    match (max_1, max_2) {
        (None, None) => None,
        // One or both reported a safety bound — take the tighter
        // (smaller) bound since both surfaces have to accommodate
        // the rolling ball simultaneously.
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
    }
}
