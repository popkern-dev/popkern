//! Tests for the region-scoped pcurve recompute and stale repair
//! (PK-274).

use crate::geometry::pcurve::Pcurve;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// PK-274: the scoped recompute over the full face set matches the
/// whole-model form (parity), and dead/duplicate handles are ignored.
#[test]
fn scoped_recompute_full_face_set_matches_whole_model() {
    let tol = Tolerance::default();
    let mut whole = BRep::unit_cube(tol).unwrap();
    let mut scoped = whole.clone();
    for he in whole.half_edges.values_mut() {
        he.pcurve = None;
    }
    for he in scoped.half_edges.values_mut() {
        he.pcurve = None;
    }
    whole.recompute_all_pcurves().unwrap();
    let mut faces: Vec<_> = scoped.faces.keys().collect();
    let dup = faces[0];
    faces.push(dup);
    scoped.recompute_pcurves_for(&faces).unwrap();
    for (heh, he) in &whole.half_edges {
        let other = &scoped.half_edges[heh];
        assert_eq!(he.pcurve.is_some(), other.pcurve.is_some());
        if let (Some(a), Some(b)) = (he.pcurve.as_ref(), other.pcurve.as_ref()) {
            let mid_a = a.point_at(0.5 * (a.t_start + a.t_end));
            let mid_b = b.point_at(0.5 * (b.t_start + b.t_end));
            assert!((mid_a.x - mid_b.x).abs() < 1e-12 && (mid_a.y - mid_b.y).abs() < 1e-12);
        }
    }
}

/// PK-274: `recompute_pcurves_where_missing` fills exactly the faces
/// with a gap and leaves complete faces' pcurves untouched.
#[test]
fn where_missing_fills_only_gapped_faces() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).unwrap();
    let gap_face = brep.faces.keys().next().unwrap();
    let loop_h = brep.faces[gap_face].outer_loop;
    let he0 = brep.loops[loop_h].half_edge;
    brep.half_edges[he0].pcurve = None;
    let untouched: Vec<_> = brep
        .half_edges
        .iter()
        .filter(|(h, he)| *h != he0 && he.pcurve.is_some())
        .map(|(h, he)| {
            let pc = he.pcurve.as_ref().unwrap();
            (h, pc.point_at(0.5 * (pc.t_start + pc.t_end)))
        })
        .collect();
    brep.recompute_pcurves_where_missing().unwrap();
    assert!(brep.half_edges[he0].pcurve.is_some(), "gap must be filled");
    for (h, mid_before) in untouched {
        let Some(face) = brep.loops.get(brep.half_edges[h].loop_).map(|l| l.face) else {
            continue;
        };
        if face == gap_face {
            continue; // same-face pcurves legitimately recomputed
        }
        let pc = brep.half_edges[h].pcurve.as_ref().unwrap();
        let mid_now = pc.point_at(0.5 * (pc.t_start + pc.t_end));
        assert!(
            (mid_now.x - mid_before.x).abs() < 1e-12 && (mid_now.y - mid_before.y).abs() < 1e-12,
            "complete faces' pcurves must be untouched"
        );
    }
}

/// PK-274: `stale_pcurve_handles` is read-only detection with the
/// same gate as `fix_stale_pcurves`; a clean cube reports none, and
/// the scoped repair clears an injected stale pcurve.
#[test]
fn stale_detection_and_scoped_repair_agree() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).unwrap();
    assert_eq!(brep.stale_pcurve_handles(), vec![]);
    let face = brep.faces.keys().next().unwrap();
    let loop_h = brep.faces[face].outer_loop;
    let he0 = brep.loops[loop_h].half_edge;
    // Inject a garbage pcurve: shifted far off the true trim.
    brep.half_edges[he0].pcurve = Some(Pcurve::linear((0.4, 0.4), (0.6, 0.6), 0.0, 1.0));
    let stale = brep.stale_pcurve_handles();
    assert_eq!(stale, vec![he0]);
    brep.fix_stale_pcurves_among(&[face]);
    assert_eq!(
        brep.stale_pcurve_handles(),
        vec![],
        "scoped repair over the owning face must clear the staleness"
    );
}
