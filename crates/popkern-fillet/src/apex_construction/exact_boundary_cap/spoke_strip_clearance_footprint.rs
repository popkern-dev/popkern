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
        let ((_, _), v_range) = exact_bounded_parameter_range(&vertices, &surface);
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
        let Some(u_range) = arc_angular_window(arc, &surface, u_periodic) else {
            continue;
        };
        strips.push(StripFootprint {
            face,
            surface,
            u_periodic,
            u_range,
            v_range,
        });
    }
    strips
}

/// The strip's angular window at its live end: project samples of the
/// end arc (exact on-surface geometry) to UV, unwrap the periodic
/// angle around the first sample, and pad the span. `None` when no
/// sample projects.
pub(super) fn arc_angular_window(
    arc: &NurbsCurve,
    surface: &SurfaceGeometry,
    u_periodic: bool,
) -> Option<(f64, f64)> {
    let mut reference = None;
    let mut lo = f64::MAX;
    let mut hi = f64::MIN;
    for k in 0..=SPOKE_SAMPLES {
        let t = k as f64 / SPOKE_SAMPLES as f64;
        let point = arc.point_at(t);
        let Some((u, _)) = project_point_to_uv(surface, &point) else {
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
    }
    reference?;
    let padding = (hi - lo) * ANGULAR_WINDOW_PADDING;
    Some((lo - padding, hi + padding))
}

/// Shift `u` by whole periods (2π) to the branch nearest `reference`.
pub(super) fn unwrap_periodic(u: f64, reference: f64) -> f64 {
    let period = std::f64::consts::TAU;
    u - ((u - reference) / period).round() * period
}

impl StripFootprint {
    /// Whether the point projects into the strip's footprint window
    /// (arc-derived angular span × vertex-bounded axial span) — the
    /// live footprint, not the infinite primitive.
    pub(super) fn contains_footprint_projection(&self, point: &Point3) -> bool {
        let Some((u, v)) = project_point_to_uv(&self.surface, point) else {
            return false;
        };
        let u = if self.u_periodic {
            unwrap_periodic(u, (self.u_range.0 + self.u_range.1) * 0.5)
        } else {
            u
        };
        u >= self.u_range.0 && u <= self.u_range.1 && v >= self.v_range.0 && v <= self.v_range.1
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
