//! Exact geometry payload for post-merge side-support trim edges.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::edge_triple::{EdgeAnchors, EdgeTripleOptions, compute_edge_geometry};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{CurveGeometry, CurveOnEdge};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

#[derive(Clone, Debug)]
pub(super) enum SideGeneratorCurveAuthoring {
    Ready(ExactSideGeneratorGeometry),
    Refused(String),
}

#[derive(Clone, Debug)]
pub(super) struct ExactSideGeneratorGeometry {
    curve: CurveOnEdge,
    cylinder_pcurve: Pcurve,
    support_pcurve: Pcurve,
}

impl SideGeneratorCurveAuthoring {
    pub(super) fn build(
        brep: &BRep,
        cylinder_face: FaceHandle,
        support_face: FaceHandle,
        start: VertexHandle,
        end: VertexHandle,
    ) -> Self {
        let Some(cylinder_surface) = brep
            .faces
            .get(cylinder_face)
            .and_then(|face| face.surface.as_ref())
            .map(|surface| &surface.geometry)
        else {
            return Self::Refused(format!("cylinder-face {cylinder_face:?} has no surface"));
        };
        let Some(support_surface) = brep
            .faces
            .get(support_face)
            .and_then(|face| face.surface.as_ref())
            .map(|surface| &surface.geometry)
        else {
            return Self::Refused(format!("support-face {support_face:?} has no surface"));
        };
        let Some(start_point) = brep.vertices.get(start).map(|vertex| vertex.position) else {
            return Self::Refused(format!("start vertex {start:?} is stale"));
        };
        let Some(end_point) = brep.vertices.get(end).map(|vertex| vertex.position) else {
            return Self::Refused(format!("end vertex {end:?} is stale"));
        };
        let result = compute_edge_geometry(
            cylinder_surface,
            support_surface,
            EdgeAnchors::Open {
                start: start_point,
                end: end_point,
            },
            None,
            &brep.tolerance,
            &EdgeTripleOptions::default(),
        );
        match result {
            Ok(triple) => {
                let start_error = triple
                    .curve_3d
                    .point_at(triple.t_start)
                    .distance_to(&start_point);
                let end_error = triple
                    .curve_3d
                    .point_at(triple.t_end)
                    .distance_to(&end_point);
                if start_error.max(end_error) > brep.tolerance.linear().max(1.0e-9) {
                    return Self::Refused(format!(
                        "edge-triple endpoint mismatch start={start_error:.3e} end={end_error:.3e}"
                    ));
                }
                let geometry = canonical_edge_curve(
                    triple.curve_3d,
                    triple.t_start,
                    triple.t_end,
                    brep.tolerance.linear(),
                );
                Self::Ready(ExactSideGeneratorGeometry {
                    curve: CurveOnEdge {
                        geometry,
                        t_start: triple.t_start,
                        t_end: triple.t_end,
                    },
                    cylinder_pcurve: triple.pcurve_a,
                    support_pcurve: triple.pcurve_b,
                })
            }
            Err(err) => Self::Refused(format!("edge-triple-refused={err}")),
        }
    }

    pub(super) fn is_ready(&self) -> bool {
        matches!(self, Self::Ready(_))
    }

    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            Self::Ready(_) => "exact_curve_authoring=ready".to_string(),
            Self::Refused(reason) => format!("exact_curve_authoring=refused reason={reason}"),
        }
    }

    pub(super) fn apply_to_support_split(
        &self,
        brep: &mut BRep,
        half_edge: HalfEdgeHandle,
        half_edge_twin: HalfEdgeHandle,
    ) -> Result<(), String> {
        let Self::Ready(geometry) = self else {
            return Err("side-generator curve authoring is not ready".to_string());
        };
        geometry.apply_to_support_split(brep, half_edge, half_edge_twin)
    }
}

impl ExactSideGeneratorGeometry {
    fn apply_to_support_split(
        &self,
        brep: &mut BRep,
        half_edge: HalfEdgeHandle,
        half_edge_twin: HalfEdgeHandle,
    ) -> Result<(), String> {
        let edge = brep
            .half_edges
            .get(half_edge)
            .ok_or_else(|| format!("missing side-generator half-edge {half_edge:?}"))?
            .edge;
        brep.edges
            .get_mut(edge)
            .ok_or_else(|| format!("missing side-generator edge {edge:?}"))?
            .curve = Some(self.curve.clone());
        brep.half_edges
            .get_mut(half_edge)
            .ok_or_else(|| format!("missing side-generator half-edge {half_edge:?}"))?
            .pcurve = Some(self.support_pcurve.clone());
        brep.half_edges
            .get_mut(half_edge_twin)
            .ok_or_else(|| format!("missing side-generator twin {half_edge_twin:?}"))?
            .pcurve = Some(self.support_pcurve.clone());
        let _planned_cylinder_pcurve = &self.cylinder_pcurve;
        Ok(())
    }
}

fn canonical_edge_curve(
    curve: CurveGeometry,
    t_start: f64,
    t_end: f64,
    tolerance: f64,
) -> CurveGeometry {
    curve
        .canonical_line_in_range((t_start, t_end), tolerance.max(1.0e-8))
        .unwrap_or(curve)
}
