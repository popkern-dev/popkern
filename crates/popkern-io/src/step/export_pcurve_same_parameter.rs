//! Export pcurves co-parameterized with their edge's 3D curve (PK-554).
//!
//! A STEP pcurve shares its surface curve's parameterization: at each
//! 3D-curve parameter, the pcurve's surface image is the curve point.
//! popkern stores exact curves and exact pcurves that need not share a
//! parameterization. A cylinder rim is a rational NURBS circle; its
//! wall trim is linear in angle. Such an edge used to be rewritten as a
//! 17-sample polyline cut from a trim (2.9e-2 sagitta on an r = 1.5
//! rim), and OpenCASCADE re-derives pcurves from that curve.
//!
//! Now the exact curve is kept. Each pcurve that is not co-parameterized
//! with it is rebuilt on the curve's own parameter
//! ([`same_parameter_repair`]), within the model's linear tolerance
//! when possible and its declared accuracy otherwise. A pcurve that
//! cannot be rebuilt keeps the direction-aligned conversion.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceOnFace};
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::export_pcurve_synthesis::step_export_pcurve_association_tolerance;
use super::same_parameter::same_parameter_repair;

/// Samples of the co-parameterization check.
const CO_PARAMETER_SAMPLES: usize = 17;

/// The NURBS pcurve to write for `pcurve` (a half-edge's trim, in
/// half-edge direction) on `surface`, parameterized like the edge's 3D
/// curve in the curve's direction.
pub(super) fn export_pcurve_on_edge_curve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceOnFace,
    pcurve: &Pcurve,
) -> NurbsCurve {
    let opposes = super::export_periodic_branches::pcurve_opposes_edge_curve(
        brep,
        half_edge,
        &surface.geometry,
        pcurve,
    );
    if let Some(rebuilt) = rebuilt_on_curve(brep, half_edge, surface, pcurve, opposes) {
        return rebuilt;
    }
    let mut nurbs = super::export_pcurve::export_pcurve_nurbs(
        pcurve,
        &surface.geometry,
        brep.tolerance.linear(),
    );
    if opposes {
        super::export_periodic_branches::reverse_pcurve_nurbs(&mut nurbs);
    }
    nurbs
}

/// The pcurve rebuilt on the edge curve's parameter, when it is not
/// already co-parameterized with it.
fn rebuilt_on_curve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceOnFace,
    pcurve: &Pcurve,
    opposes: bool,
) -> Option<NurbsCurve> {
    let edge = brep.half_edges.get(half_edge)?.edge;
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let linear = brep.tolerance.linear();
    if co_parameter_gap(pcurve, surface, curve, opposes) <= linear {
        return None;
    }
    let repair = same_parameter_repair(pcurve, curve, surface, true, linear).or_else(|| {
        same_parameter_repair(
            pcurve,
            curve,
            surface,
            true,
            step_export_pcurve_association_tolerance(brep),
        )
    })?;
    match repair.pcurve.geometry {
        PcurveGeometry::Nurbs(fit) => Some(fit),
        _ => None,
    }
}

/// Largest distance between the pcurve's surface image and the edge
/// curve at equal fractions, the pcurve taken in the curve's direction.
fn co_parameter_gap(
    pcurve: &Pcurve,
    surface: &SurfaceOnFace,
    curve: &CurveOnEdge,
    opposes: bool,
) -> f64 {
    let last = (CO_PARAMETER_SAMPLES - 1) as f64;
    (0..CO_PARAMETER_SAMPLES)
        .map(|k| {
            let f = k as f64 / last;
            let pf = if opposes { 1.0 - f } else { f };
            let uv = pcurve.point_at(pcurve.t_start + pf * (pcurve.t_end - pcurve.t_start));
            let on_surface = surface.geometry.point_at(uv.x, uv.y);
            let on_curve = curve
                .geometry
                .point_at(curve.t_start + f * (curve.t_end - curve.t_start));
            on_surface.distance_to(&on_curve)
        })
        .fold(0.0_f64, f64::max)
}
