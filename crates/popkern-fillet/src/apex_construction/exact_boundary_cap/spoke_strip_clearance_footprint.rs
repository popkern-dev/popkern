//! Strip footprint-window construction for the spoke clearance pass.
//!
//! Nested inside `spoke_strip_clearance` (shares its private types);
//! extracted to hold the parent under the 500-line module budget.

use popkern_core::geometry::analytic_surfaces::analytic_signed_distance;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::{exact_bounded_parameter_range, project_point_to_uv};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::BRep;

use super::super::corner::PreparedExactBoundarySectorBoundary;
use super::super::sector::ExactBoundaryCapSectorInput;
use super::{SPOKE_SAMPLES, StripFootprint};

/// Fractional padding applied to the arc-derived angular window, so
/// the clearance zone covers the strip's rails plus a small margin.
const ANGULAR_WINDOW_PADDING: f64 = 0.05;

/// Gather the corner's unique fillet-strip supports with their
/// footprint windows. Strips whose surface has no implicit form or
/// whose end arc cannot be projected are skipped — the pass corrects
/// only what it can measure honestly.
pub(super) fn collect_strip_footprints(
    brep: &BRep,
    sectors: &[ExactBoundaryCapSectorInput],
    prepared: &[PreparedExactBoundarySectorBoundary],
) -> Vec<StripFootprint> {
    let mut strips: Vec<StripFootprint> = Vec::new();
    for sector in sectors {
        let face = sector.fillet_support_face;
        if strips.iter().any(|strip| strip.face == face) {
            continue;
        }
        let Some(surface_on_face) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let surface = surface_on_face.geometry.clone();
        let origin_probe = surface.point_at(0.0, 0.0);
        if analytic_signed_distance(&surface, &origin_probe).is_none() {
            continue;
        }
        let vertices = brep.face_outer_loop_positions(face);
        if vertices.is_empty() {
            continue;
        }
        let ((_, _), vertex_v) = exact_bounded_parameter_range(&vertices, &surface);
        let Some(arc) = prepared
            .iter()
            .find(|prep| prep.sector_index == sector.sector_index)
            .map(|prep| &prep.fillet_boundary.curve)
        else {
            continue;
        };
        let u_periodic = matches!(
            surface,
            SurfaceGeometry::Cylinder { .. }
                | SurfaceGeometry::Cone { .. }
                | SurfaceGeometry::Sphere { .. }
                | SurfaceGeometry::Torus { .. }
        );
        let Some((u_range, arc_v)) = arc_windows(arc, &surface, u_periodic) else {
            continue;
        };
        let Some(arc_profile) = arc_uv_profile(arc, &surface, u_periodic) else {
            continue;
        };
        let arc_mid = (arc_v.0 + arc_v.1) * 0.5;
        let far_side_positive = (vertex_v.1 - arc_mid) >= (arc_mid - vertex_v.0);
        strips.push(StripFootprint {
            face,
            surface,
            u_periodic,
            u_range,
            v_range: corner_anchored_axial_window(arc_v, vertex_v),
            arc_profile,
            far_side_positive,
        });
    }
    strips
}

/// The end arc as a `(unwrapped u, v)` profile sorted by u — the
/// footprint's true corner-side boundary curve for
/// [`StripFootprint::corner_side_of_arc`].
fn arc_uv_profile(
    arc: &NurbsCurve,
    surface: &SurfaceGeometry,
    u_periodic: bool,
) -> Option<Vec<(f64, f64)>> {
    let mut profile: Vec<(f64, f64)> = Vec::with_capacity(SPOKE_SAMPLES + 1);
    let mut reference = None;
    for k in 0..=SPOKE_SAMPLES {
        let t = k as f64 / SPOKE_SAMPLES as f64;
        let point = arc.point_at(t);
        let Some((u, v)) = project_point_to_uv(surface, &point) else {
            continue;
        };
        let u = match reference {
            None => {
                reference = Some(u);
                u
            }
            Some(reference) if u_periodic => unwrap_periodic(u, reference),
            Some(_) => u,
        };
        profile.push((u, v));
    }
    if profile.len() < 2 {
        return None;
    }
    profile.sort_by(|a, b| a.0.total_cmp(&b.0));
    Some(profile)
}

/// The strip's windows at its live end: project samples of the end
/// arc (exact on-surface geometry) to UV, unwrap the periodic angle
/// around the first sample, and pad the angular span. Returns the
/// padded angular window and the arc's raw axial (v) span. `None`
/// when no sample projects.
pub(super) fn arc_windows(
    arc: &NurbsCurve,
    surface: &SurfaceGeometry,
    u_periodic: bool,
) -> Option<((f64, f64), (f64, f64))> {
    let mut reference = None;
    let mut lo = f64::MAX;
    let mut hi = f64::MIN;
    let mut v_lo = f64::MAX;
    let mut v_hi = f64::MIN;
    for k in 0..=SPOKE_SAMPLES {
        let t = k as f64 / SPOKE_SAMPLES as f64;
        let point = arc.point_at(t);
        let Some((u, v)) = project_point_to_uv(surface, &point) else {
            continue;
        };
        let u = match reference {
            None => {
                reference = Some(u);
                u
            }
            Some(reference) if u_periodic => unwrap_periodic(u, reference),
            Some(_) => u,
        };
        lo = lo.min(u);
        hi = hi.max(u);
        v_lo = v_lo.min(v);
        v_hi = v_hi.max(v);
    }
    reference?;
    let padding = (hi - lo) * ANGULAR_WINDOW_PADDING;
    Some(((lo - padding, hi + padding), (v_lo, v_hi)))
}

/// Fraction of the strip's axial extent (beyond the corner arc) the
/// clearance zone covers. The FAR half of the strip is untrustworthy
/// mid-candidate — its true trim is the oblique end section that
/// `strip_end_arcs` authors only in the post-pipeline, so at spoke
/// time the far boundary is still the Phase C incoherent state
/// (measured pyramid4_mixed: far-side "dips" at v beyond the real far
/// trim, blessed by the honest validator, regressed hover/repair pins
/// when corrected). The corner band is where the class this pass
/// exists for lives: pyramid3_extreme's real dips hug the corner arc
/// at +0.0019 axially.
const CORNER_BAND_FRACTION: f64 = 0.25;

/// Anchor the axial window's CORNER end to the end arc and extend it
/// only a corner-band fraction of the way toward the far end. The
/// vertex bound alone over-includes pole-ward of the arc (off-surface
/// loop vertices) and reaches into the untrustworthy far half; the
/// corner side is whichever end of the vertex window the arc sits
/// nearest, snapped to the arc's own span.
fn corner_anchored_axial_window(arc_v: (f64, f64), vertex_v: (f64, f64)) -> (f64, f64) {
    let arc_mid = (arc_v.0 + arc_v.1) * 0.5;
    if (arc_mid - vertex_v.0).abs() <= (arc_mid - vertex_v.1).abs() {
        let far = vertex_v.1.max(arc_v.1);
        (arc_v.0, arc_v.1 + (far - arc_v.1) * CORNER_BAND_FRACTION)
    } else {
        let far = vertex_v.0.min(arc_v.0);
        (arc_v.0 - (arc_v.0 - far) * CORNER_BAND_FRACTION, arc_v.1)
    }
}

/// Shift `u` by whole periods (2π) to the branch nearest `reference`.
pub(super) fn unwrap_periodic(u: f64, reference: f64) -> f64 {
    let period = std::f64::consts::TAU;
    u - ((u - reference) / period).round() * period
}

impl StripFootprint {
    /// Whether the point projects into the strip's footprint window:
    /// arc-derived angular span, vertex-bounded axial span, AND on the
    /// strip's far side of the end-arc profile at the point's own
    /// angular position — the live footprint, not the infinite
    /// primitive and not a rectangle over an oblique end arc.
    pub(super) fn contains_footprint_projection(&self, point: &Point3) -> bool {
        let Some((u, v)) = project_point_to_uv(&self.surface, point) else {
            return false;
        };
        let u = if self.u_periodic {
            unwrap_periodic(u, (self.u_range.0 + self.u_range.1) * 0.5)
        } else {
            u
        };
        u >= self.u_range.0
            && u <= self.u_range.1
            && v >= self.v_range.0
            && v <= self.v_range.1
            && self.corner_side_of_arc(u, v)
    }

    /// Whether `(u, v)` lies on the strip's interior side of the end
    /// arc at angular position `u` (linear interpolation over the arc
    /// profile, clamped at its ends).
    fn corner_side_of_arc(&self, u: f64, v: f64) -> bool {
        let profile = &self.arc_profile;
        let arc_v = if u <= profile[0].0 {
            profile[0].1
        } else if u >= profile[profile.len() - 1].0 {
            profile[profile.len() - 1].1
        } else {
            let mut value = profile[0].1;
            for pair in profile.windows(2) {
                let (u0, v0) = pair[0];
                let (u1, v1) = pair[1];
                if u >= u0 && u <= u1 {
                    let alpha = if u1 > u0 { (u - u0) / (u1 - u0) } else { 0.0 };
                    value = v0 + (v1 - v0) * alpha;
                    break;
                }
            }
            value
        };
        if self.far_side_positive {
            v >= arc_v
        } else {
            v <= arc_v
        }
    }

    /// Whether the point lies within the strip's axial (v) window,
    /// regardless of angular position — the zone where the primitive's
    /// material side is meaningful for the rim side vote.
    pub(super) fn contains_axial_projection(&self, point: &Point3) -> bool {
        let Some((_, v)) = project_point_to_uv(&self.surface, point) else {
            return false;
        };
        v >= self.v_range.0 && v <= self.v_range.1
    }
}
