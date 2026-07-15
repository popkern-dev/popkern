//! Trace-gated per-seam provenance census for cap/support seams —
//! how far each live seam actor (endpoint vertices, edge curve,
//! support-surface trace) sits from the cap surface.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};
use std::collections::BTreeSet;

/// Trace-gated per-seam provenance census: how far each live seam
/// actor (endpoint vertices, edge curve, support-surface trace) sits
/// from the cap surface — separates stale-actor classes before any
/// repair mutates them.
pub(super) fn trace_cap_support_seam_census(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_faces: &BTreeSet<FaceHandle>,
) {
    if !crate::direct_apex_trace::enabled() {
        return;
    }
    let Some(cap_surface) = brep
        .faces
        .get(cap_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return;
    };
    let cap_dist = |p: &Point3| -> f64 {
        let uv = super::cap_uv_point(&cap_surface, p);
        cap_surface.point_at(uv.x, uv.y).distance_to(p)
    };
    for half_edge in super::cap_face_half_edges(brep, cap_face) {
        let Some(support_face) = super::twin_face_for_half_edge(brep, half_edge) else {
            continue;
        };
        if cap_faces.contains(&support_face) {
            continue;
        }
        let Some(he) = brep.half_edges.get(half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(v0), Some(v1)) = (brep.vertices.get(he.origin), brep.vertices.get(twin.origin))
        else {
            continue;
        };
        let support_surface = brep
            .faces
            .get(support_face)
            .and_then(|face| face.surface.as_ref())
            .map(|surface| surface.geometry.clone());
        let curve = brep.edges.get(he.edge).and_then(|edge| edge.curve.as_ref());
        let curve_mid_caps: Vec<f64> = curve
            .map(|c| {
                [0.25, 0.5, 0.75]
                    .iter()
                    .map(|&f| {
                        let t = c.t_start + (c.t_end - c.t_start) * f;
                        let p = popkern_core::geometry::curve::CurveEval::point_at(&c.geometry, t);
                        cap_dist(&p)
                    })
                    .collect()
            })
            .unwrap_or_default();
        let curve_mid_support: Option<f64> = curve.zip(support_surface.as_ref()).map(|(c, ss)| {
            let t = c.t_start + (c.t_end - c.t_start) * 0.5;
            let p = popkern_core::geometry::curve::CurveEval::point_at(&c.geometry, t);
            let uv = super::unclamped_uv_point(ss, &p);
            ss.point_at(uv.x, uv.y).distance_to(&p)
        });
        eprintln!(
            "[seam-census] he {half_edge:?} edge {:?} support {support_face:?}: v0_cap={:.3e} v1_cap={:.3e} curve_cap={:?} curve_support={:?}",
            he.edge,
            cap_dist(&v0.position),
            cap_dist(&v1.position),
            curve_mid_caps
                .iter()
                .map(|d| format!("{d:.3e}"))
                .collect::<Vec<_>>(),
            curve_mid_support.map(|d| format!("{d:.3e}")),
        );
    }
}
