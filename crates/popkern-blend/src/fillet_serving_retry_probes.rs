//! The served-body probes the serving retry acts on (PK-250, PK-260,
//! PK-597): per corner, whether the first pass's direct-apex serve
//! self-intersects, declined at the serving-shape gate, or creases
//! across its own cap seams.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::blend_continuity_evidence::{ContinuityGates, measure_blend_continuity_evidence};
use crate::types::{BlendResult, ContinuityLevel};

/// Samples per seam for the served crease check (the evidence default).
const CREASE_SAMPLES_PER_SEAM: usize = 33;

/// Corners whose direct-apex construction declined at the
/// serving-shape gate (PK-260), with the recorded decline reason. The
/// mid-pipeline fallback serves a degraded vertex-blend body (direct-apex
/// anchor convention); the exclusion re-run serves clean.
pub(super) fn serving_shape_declined_corners(
    result: &BlendResult,
) -> Vec<(popkern_core::topology::VertexHandle, String)> {
    result
        .repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            route.route == crate::types::ApexConstructionRoute::DirectApexFallbackToVertexBlend
        })
        .filter_map(|route| {
            route
                .fallback_reason
                .as_ref()
                .filter(|reason| reason.contains(crate::types::SERVING_SHAPE_DECLINE_PREFIX))
                .map(|reason| (route.vertex, format!("serving-shape decline: {reason}")))
        })
        .collect()
}

/// Served direct-apex corners whose region has a face whose loops
/// bound no region of its periodic chart (an edge run the long way
/// round, PK-792: a cube corner cap's rim wrapped three quarters round
/// its r = 0.1 strip). Structural, so the corner takes the exclusion
/// re-run directly: a footprint repair moves cap interiors only and
/// cannot shorten a rim.
pub(super) fn served_direct_apex_unbounded_trims(
    result: &BlendResult,
) -> Vec<(popkern_core::topology::VertexHandle, String)> {
    result
        .repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            matches!(route.route, crate::types::ApexConstructionRoute::DirectApex)
                && !route.built_faces.is_empty()
        })
        .filter_map(|route| {
            let region = crate::apex_construction::exact_boundary_cap::cap_quality_postcondition::
                candidate_region(&result.brep, &route.built_faces);
            result
                .brep
                .validate_loop_period_winding_among(&region)
                .first()
                .map(|error| {
                    (
                        route.vertex,
                        format!("served trim bounds no region: {error}"),
                    )
                })
        })
        .collect()
}

/// Served direct-apex corners whose region self-intersects, with the
/// validator's face-pair detail.
pub(super) fn served_direct_apex_self_intersections(
    result: &BlendResult,
) -> Vec<(popkern_core::topology::VertexHandle, String)> {
    result
        .repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            matches!(route.route, crate::types::ApexConstructionRoute::DirectApex)
                && !route.built_faces.is_empty()
        })
        .filter_map(|route| {
            crate::apex_construction::exact_boundary_cap::cap_quality_postcondition::
                validate_candidate_region_self_intersection(&result.brep, &route.built_faces)
                .err()
                .map(|detail| (route.vertex, detail))
        })
        .collect()
}

/// Served direct-apex corners whose cap creases across its own seams
/// (PK-597), with the worst seam. Measured by the whole-seam continuity
/// evidence on the served body, restricted to seams between two of the
/// corner's cap faces — every such seam, whatever its angle: a cap
/// has no designed internal crease (PK-609; under the blend-to-blend
/// angle bound the cube 0.05/0.2/0.3 cap's sector seam, certified
/// normal sine 0.9998, read as a designed crease and never fired).
/// The PK-197 mid-seam instrument is not used: it reads fractions
/// 0.4-0.6 only (1.5e-8 on a mixed n = 4 G1 apex whose cap creases
/// 28.6° toward the seam ends, PK-607). A cap-cap normal jump above the
/// commercial G1 gate breaks the G1 every corner blend claims, so the
/// corner takes the exclusion re-run (measured: the n = 4 apex then
/// serves the setback rung, certified G1 at 2.8e-9).
pub(super) fn served_direct_apex_creases(
    result: &BlendResult,
) -> Vec<(popkern_core::topology::VertexHandle, String)> {
    let gates = ContinuityGates::commercial();
    result
        .repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            matches!(route.route, crate::types::ApexConstructionRoute::DirectApex)
                && route.built_faces.len() > 1
        })
        .filter_map(|route| {
            let evidence = measure_blend_continuity_evidence(
                &result.brep,
                &route.built_faces,
                ContinuityLevel::G1,
                &gates,
                CREASE_SAMPLES_PER_SEAM,
            );
            let worst = evidence
                .seams
                .iter()
                .filter(|seam| both_sides_in(&result.brep, seam.edge, &route.built_faces))
                .map(|seam| seam.worst_g1_normal)
                .fold(0.0_f64, f64::max);
            (worst > gates.g1_normal).then(|| {
                (
                    route.vertex,
                    format!(
                        "served cap-cap crease: normal sine {worst:.3e} over the G1 gate {:e}",
                        gates.g1_normal
                    ),
                )
            })
        })
        .collect()
}

/// Whether both faces across `edge` belong to `faces`.
fn both_sides_in(brep: &BRep, edge: EdgeHandle, faces: &[FaceHandle]) -> bool {
    let sides: Vec<FaceHandle> = brep
        .half_edges
        .values()
        .filter(|he| he.edge == edge)
        .filter_map(|he| brep.loops.get(he.loop_).map(|l| l.face))
        .collect();
    sides.len() == 2 && sides.iter().all(|f| faces.contains(f))
}
