use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{
    BRep, BrepEdgeSurfaceCompatibilityOptions, analyze_brep_edge_surface_compatibility,
};

pub(crate) fn pcurve_side_agreement_summary(brep: &BRep, limit: usize) -> Vec<String> {
    let edge_surface = analyze_brep_edge_surface_compatibility(
        brep,
        BrepEdgeSurfaceCompatibilityOptions {
            tolerance: 1e-6,
            max_anchor_motion: 0.0,
            samples: 9,
            require_all_compatible: false,
            validate_output: false,
        },
    );
    let mut out = Vec::new();
    for finding in edge_surface.findings.iter().take(limit) {
        let edge_handle = finding.edge;
        let Some(edge) = brep.edges.get(edge_handle) else {
            continue;
        };
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let Some(first_he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let pair = [edge.half_edge, first_he.twin];
        let mut samples = Vec::new();
        for heh in pair {
            let Some(he) = brep.half_edges.get(heh) else {
                continue;
            };
            let Some(pcurve) = he.pcurve.as_ref() else {
                continue;
            };
            let Some(loop_) = brep.loops.get(he.loop_) else {
                continue;
            };
            let Some(surface) = brep
                .faces
                .get(loop_.face)
                .and_then(|face| face.surface.as_ref())
            else {
                continue;
            };
            let t = (pcurve.t_start + pcurve.t_end) * 0.5;
            let uv = pcurve.point_at(t);
            samples.push((loop_.face, surface.geometry.point_at(uv.x, uv.y)));
        }
        if samples.len() == 2 {
            let t_mid = (curve.t_start + curve.t_end) * 0.5;
            let p_edge = curve.geometry.point_at(t_mid);
            out.push(format!(
                "{edge_handle:?}: pcurve_pair_gap={:.6} edge_gap_a={:.6} \
                 edge_gap_b={:.6} {:?}/{:?}",
                samples[0].1.distance_to(&samples[1].1),
                p_edge.distance_to(&samples[0].1),
                p_edge.distance_to(&samples[1].1),
                samples[0].0,
                samples[1].0
            ));
        }
    }
    out
}
