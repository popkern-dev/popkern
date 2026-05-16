use super::*;
use crate::geometry::{NurbsCurve, Point3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn options() -> BrepTrimRepairOptions {
    BrepTrimRepairOptions {
        endpoint_tolerance: 1e-6,
        sampled_tolerance: 1e-6,
        samples: 5,
        validate_output: true,
    }
}

#[test]
fn trim_diagnostics_accept_clean_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let report = analyze_brep_trims(&brep, options());

    assert!(report.is_clean());
    assert!(report.inspected_half_edges > 0);
    assert_eq!(report.missing_pcurves, 0);
    assert_eq!(report.stale_endpoint_pcurves, 0);
    assert_eq!(report.sampled_mismatches, 0);
}

#[test]
fn migration_matrix_trim_report_contract() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let half_edges = face_half_edges_for_test(&brep, face);

    let clean = analyze_brep_trims(&brep, options());
    assert_eq!(clean.missing_pcurves, 0);
    assert_eq!(clean.stale_endpoint_pcurves, 0);
    assert_eq!(clean.sampled_mismatches, 0);

    brep.half_edges[half_edges[0]].pcurve = None;
    let missing = analyze_brep_trims(&brep, options());
    assert_eq!(missing.missing_pcurves, 1);

    corrupt_pcurve(&mut brep, half_edges[1]);
    let stale = analyze_brep_trims(&brep, options());
    assert!(stale.stale_endpoint_pcurves > 0);
    assert!(stale.max_endpoint_deviation > 1.0);
}

#[test]
fn trim_repair_fills_missing_pcurves_transactionally() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let half_edges = face_half_edges_for_test(&brep, face);
    for half_edge in &half_edges {
        brep.half_edges[*half_edge].pcurve = None;
    }

    let before = analyze_brep_trims(&brep, options());
    assert_eq!(before.missing_pcurves, half_edges.len());

    let repaired = repair_brep_trims_transactionally(&brep, options()).unwrap();

    assert!(brep.half_edges[half_edges[0]].pcurve.is_none());
    assert!(repaired.report.after.is_clean());
    assert_eq!(repaired.report.repaired_faces, 1);
    assert_eq!(repaired.report.repaired_half_edges, half_edges.len());
    for half_edge in half_edges {
        assert!(repaired.brep.half_edges[half_edge].pcurve.is_some());
    }
}

#[test]
fn trim_repair_replaces_corrupted_pcurves_transactionally() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let half_edges = face_half_edges_for_test(&brep, face);
    for half_edge in &half_edges {
        corrupt_pcurve(&mut brep, *half_edge);
    }

    let before = analyze_brep_trims(&brep, options());
    assert!(before.stale_endpoint_pcurves > 0);
    assert!(before.max_endpoint_deviation > 1.0);

    let repaired = repair_brep_trims_transactionally(&brep, options()).unwrap();

    assert!(
        analyze_brep_trims(&brep, options()).max_endpoint_deviation > 1.0,
        "input BRep remains corrupted"
    );
    assert!(repaired.report.after.is_clean());
    assert!(repaired.report.after.max_endpoint_deviation < 1e-6);
}

#[test]
fn trim_repair_noops_clean_cube_without_mutating_counts() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let repaired = repair_brep_trims_transactionally(&brep, options()).unwrap();

    assert_eq!(repaired.report.repaired_faces, 0);
    assert_eq!(repaired.report.repaired_half_edges, 0);
    assert!(repaired.report.before.is_clean());
    assert!(repaired.report.after.is_clean());
    assert_eq!(repaired.brep.faces.len(), brep.faces.len());
    assert_eq!(repaired.brep.half_edges.len(), brep.half_edges.len());
}

fn face_half_edges_for_test(
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

fn corrupt_pcurve(brep: &mut BRep, half_edge: crate::topology::HalfEdgeHandle) {
    let Some(pcurve) = brep.half_edges[half_edge].pcurve.as_ref() else {
        return;
    };
    let pcurve = pcurve.to_nurbs();
    let bad = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(pcurve.control_points[0].x + 100.0, 100.0, 0.0),
            Point3::new(100.0, pcurve.control_points[1].y + 100.0, 0.0),
        ],
        weights: vec![1.0, 1.0],
    };
    brep.half_edges[half_edge].pcurve = Some(bad.into());
}
