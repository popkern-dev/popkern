//! Non-Gregory fallback sampling for the KP20 cap/support boundary solve.
//!
//! When the cap face surface is not `RectangularBiNonic` (PnS produces
//! `SurfaceGeometry::Nurbs`) or the cap pcurve is routed off the
//! canonical UV boundary, the boundary-row fit in the main solve
//! module cannot run — both the row-fit step and the canonical
//! pcurve span check are `RectangularBiNonic`-specific by design.
//! The canonical *edge curve* reauthoring step is not: it just needs
//! a half-edge and support-derived 3D target samples to drive
//! [`crate::construct::kp20_polar_cap_support_edge_curve::reauthor_from_samples`].
//!
//! This module isolates that fallback so the main solver file stays
//! under the 500-line structural cap.

use popkern_core::geometry::Point3;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::BoundarySide;
use super::geometry;

/// Number of samples per edge — kept in sync with `SAMPLE_COUNT`
/// in the parent module. Both paths use the same density so the
/// resulting canonical edge curves have the same support-side
/// fidelity regardless of whether the row-fit ran.
pub(super) const SAMPLE_COUNT: usize = 12;

/// Output shape consumed by the parent module's `reauthor_from_samples`
/// loop. The non-row-fit fields (`cap_side`, `cap_span_start`,
/// `cap_span_end`, `cap_forward`) are placeholders — the parent's
/// `external_cap_support_edge_gap` helper short-circuits on the
/// absent `RectangularBiNonic` patch, so they are never read on
/// this path.
#[derive(Debug, Clone)]
pub(super) struct EdgeCurveOnlySamples {
    pub(super) cap_half_edge: HalfEdgeHandle,
    pub(super) cap_side: BoundarySide,
    pub(super) cap_span_start: f64,
    pub(super) cap_span_end: f64,
    pub(super) cap_forward: bool,
    pub(super) targets: Vec<Point3>,
}

/// Mirrors the support-side sampling that the Gregory row-fit path
/// performs, but stops short of populating the row-fit data
/// structures. Endpoints are pinned to the half-edge's origin
/// vertex positions (the topology authority); interior samples
/// come from `support_surface.point_at(support_pcurve(t))`.
pub(super) fn collect_edge_curve_samples_only(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    support_surface: &SurfaceGeometry,
    support_pcurve: &Pcurve,
) -> Option<EdgeCurveOnlySamples> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let twin_data = brep.half_edges.get(twin)?;
    let cap_start = brep.vertices[half_edge_data.origin].position;
    let cap_end = brep.vertices[twin_data.origin].position;
    let support_forward =
        geometry::support_pcurve_forward(support_surface, support_pcurve, cap_start, cap_end);
    let (support_t0, support_t1) = support_pcurve.parameter_range();
    let mut targets = Vec::with_capacity(SAMPLE_COUNT);
    for index in 0..SAMPLE_COUNT {
        let fraction = index as f64 / (SAMPLE_COUNT - 1) as f64;
        let support_fraction = if support_forward {
            fraction
        } else {
            1.0 - fraction
        };
        let support_parameter = geometry::lerp_scalar(support_t0, support_t1, support_fraction);
        let target = if index == 0 {
            cap_start
        } else if index + 1 == SAMPLE_COUNT {
            cap_end
        } else {
            geometry::support_surface_point(support_surface, support_pcurve, support_parameter)
        };
        targets.push(target);
    }
    Some(EdgeCurveOnlySamples {
        cap_half_edge: half_edge,
        cap_side: BoundarySide::UMin,
        cap_span_start: 0.0,
        cap_span_end: 1.0,
        cap_forward: true,
        targets,
    })
}
