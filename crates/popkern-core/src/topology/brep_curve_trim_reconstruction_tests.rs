use super::*;
use crate::geometry::trim_domain::TrimDomainCandidateKind;
use crate::geometry::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn options(replace_existing: bool) -> BrepCurveTrimReconstructionOptions {
    BrepCurveTrimReconstructionOptions {
        tolerance: 1.0e-6,
        samples: 9,
        replace_existing,
        replace_stale_existing: false,
        expand_bounded_surface_domains: false,
        require_all_reconstructed: false,
        validate_output: true,
    }
}

fn stale_options() -> BrepCurveTrimReconstructionOptions {
    BrepCurveTrimReconstructionOptions {
        replace_stale_existing: true,
        ..options(false)
    }
}

#[test]
fn reconstructs_missing_trim_from_existing_edge_curve_transactionally() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let half_edge = surfaced_half_edge(&brep);
    bound_owning_face_to_existing_uv_box(&mut brep, half_edge);
    brep.half_edges[half_edge].pcurve = None;

    let result = reconstruct_curve_backed_trims_transactionally(&brep, options(false)).unwrap();

    assert!(brep.half_edges[half_edge].pcurve.is_none());
    assert!(result.brep.half_edges[half_edge].pcurve.is_some());
    assert_eq!(result.report.targeted_half_edges, 1);
    assert_eq!(result.report.reconstructed_half_edges, 1);
    assert_eq!(result.report.declined_half_edges, 0);
    assert!(
        !result.report.shared_fact_ids.is_empty(),
        "curve-trim reconstruction report must expose source shared fact ids"
    );
    assert!(matches!(
        result
            .report
            .decisions
            .iter()
            .find(|d| d.half_edge == half_edge),
        Some(BrepCurveTrimReconstructionDecision {
            outcome: BrepCurveTrimReconstructionOutcome::Reconstructed {
                kind: TrimDomainCandidateKind::BoundarySubspan { .. }
            },
            ..
        })
    ));
    let decision = result
        .report
        .decisions
        .iter()
        .find(|d| d.half_edge == half_edge)
        .expect("target decision exists");
    assert!(
        !decision.shared_fact_ids.is_empty(),
        "reconstructed trim decision must reference shared facts"
    );
}

#[test]
fn replaces_corrupted_existing_trim_when_requested() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let half_edge = surfaced_half_edge(&brep);
    bound_owning_face_to_existing_uv_box(&mut brep, half_edge);
    brep.half_edges[half_edge].pcurve = Some(crate::geometry::pcurve::Pcurve::linear(
        (0.2, 0.2),
        (0.8, 0.8),
        0.0,
        1.0,
    ));

    let result = reconstruct_curve_backed_trims_transactionally(&brep, options(true)).unwrap();
    let pcurve = result.brep.half_edges[half_edge].pcurve.as_ref().unwrap();
    let surface = &result.brep.faces
        [result.brep.loops[result.brep.half_edges[half_edge].loop_].face]
        .surface
        .as_ref()
        .unwrap()
        .geometry;
    let start_uv = pcurve.point_at(pcurve.t_start);
    let start = surface.point_at(start_uv.x, start_uv.y);

    assert!(
        start.distance_to(&result.brep.vertices[result.brep.half_edges[half_edge].origin].position)
            < 1.0e-6
    );
    assert!(result.report.reconstructed_half_edges > 0);
}

#[test]
fn replaces_only_stale_existing_trim_when_enabled() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let half_edge = surfaced_half_edge(&brep);
    bound_owning_face_to_existing_uv_box(&mut brep, half_edge);
    let clean_count = brep
        .half_edges
        .values()
        .filter(|he| he.pcurve.is_some())
        .count();
    brep.half_edges[half_edge].pcurve = Some(crate::geometry::pcurve::Pcurve::linear(
        (0.2, 0.2),
        (0.8, 0.8),
        0.0,
        1.0,
    ));

    let result = reconstruct_curve_backed_trims_transactionally(&brep, stale_options()).unwrap();

    assert_eq!(result.report.targeted_half_edges, 1);
    assert_eq!(result.report.reconstructed_half_edges, 1);
    assert_eq!(
        result
            .report
            .decisions
            .iter()
            .filter(|decision| matches!(
                decision.outcome,
                BrepCurveTrimReconstructionOutcome::RetainedExisting
            ))
            .count(),
        clean_count - 1
    );
}

#[test]
fn declines_missing_edge_curve_without_mutating_existing_clone() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let half_edge = surfaced_half_edge(&brep);
    let edge = brep.half_edges[half_edge].edge;
    bound_owning_face_to_existing_uv_box(&mut brep, half_edge);
    brep.half_edges[half_edge].pcurve = None;
    brep.edges[edge].curve = None;

    let result = reconstruct_curve_backed_trims_transactionally(&brep, options(false)).unwrap();
    let decision = result
        .report
        .decisions
        .iter()
        .find(|decision| decision.half_edge == half_edge)
        .unwrap();

    assert!(result.brep.half_edges[half_edge].pcurve.is_none());
    assert!(
        !decision.shared_fact_ids.is_empty(),
        "declined reconstruction decision must reference shared facts"
    );
    assert_eq!(
        decision.outcome,
        BrepCurveTrimReconstructionOutcome::Declined(
            BrepCurveTrimReconstructionDeclineReason::MissingEdgeCurve
        )
    );
}

fn surfaced_half_edge(brep: &BRep) -> crate::topology::HalfEdgeHandle {
    brep.half_edges
        .iter()
        .find_map(|(half_edge, data)| {
            let face = brep.faces.get(brep.loops.get(data.loop_)?.face)?;
            (face.surface.is_some() && brep.edges[data.edge].curve.is_some()).then_some(half_edge)
        })
        .unwrap()
}

fn bound_owning_face_to_existing_uv_box(
    brep: &mut BRep,
    half_edge: crate::topology::HalfEdgeHandle,
) {
    let face = brep.loops[brep.half_edges[half_edge].loop_].face;
    let half_edges = face_half_edges(brep, face);
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;
    for he in half_edges {
        let Some(pcurve) = brep.half_edges[he].pcurve.as_ref() else {
            continue;
        };
        for t in [pcurve.t_start, pcurve.t_end] {
            let uv = pcurve.point_at(t);
            u_min = u_min.min(uv.x);
            u_max = u_max.max(uv.x);
            v_min = v_min.min(uv.y);
            v_max = v_max.max(uv.y);
        }
    }
    let surface = brep.faces[face].surface.as_ref().unwrap().geometry.clone();
    brep.faces[face].surface.as_mut().unwrap().geometry =
        SurfaceGeometry::Nurbs(surface.to_nurbs_in_range((u_min, u_max), (v_min, v_max)));
}

fn face_half_edges(
    brep: &BRep,
    face: crate::topology::FaceHandle,
) -> Vec<crate::topology::HalfEdgeHandle> {
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    loop {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
    }
    out
}
