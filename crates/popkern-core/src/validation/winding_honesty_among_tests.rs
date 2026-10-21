//! The region audit gives the whole-body audit's verdicts (PK-815).

use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};

fn bodies() -> Vec<(&'static str, BRep)> {
    let tol = Tolerance::default();
    vec![
        ("cube", BRep::unit_cube(tol).unwrap()),
        ("cylinder", BRep::cylinder(1.0, 2.0, tol).unwrap()),
        ("cone", BRep::cone(1.0, 2.0, tol).unwrap()),
        ("sphere", BRep::sphere(1.0, tol).unwrap()),
        ("torus", BRep::torus(3.0, 1.0, tol).unwrap()),
    ]
}

/// A torus whose faces all lost their outer-loop pcurves but `anchor`:
/// every other verdict comes from a tangency-witness chain rooted there
/// (the grid's faces meet G1, so every neighbour can vote).
fn chained_torus(anchor: usize, flip: Option<usize>) -> BRep {
    let mut b = BRep::torus(3.0, 1.0, Tolerance::default()).unwrap();
    let faces: Vec<FaceHandle> = b.faces.keys().collect();
    for (k, &f) in faces.iter().enumerate() {
        if Some(k) == flip
            && let Some(s) = b.faces[f].surface.as_mut()
        {
            s.orientation = !s.orientation;
        }
        if k == anchor {
            continue;
        }
        let start = b.loops[b.faces[f].outer_loop].half_edge;
        let mut cur = start;
        loop {
            b.half_edges[cur].pcurve = None;
            cur = b.half_edges[cur].next;
            if cur == start {
                break;
            }
        }
    }
    b
}

fn assert_region_matches_whole(label: &str, b: &BRep) {
    let whole = b.audit_winding_honesty();
    for v in &whole {
        let among = b.audit_winding_honesty_among(&[v.face]);
        assert_eq!(
            among.len(),
            1,
            "{label}: face {:?} missing from its own audit",
            v.face
        );
        assert_eq!(
            (among[0].honest, among[0].reason),
            (v.honest, v.reason),
            "{label}: face {:?} region verdict differs",
            v.face
        );
    }
    let all: Vec<FaceHandle> = b.faces.keys().collect();
    let among = b.audit_winding_honesty_among(&all);
    let key = |vs: &[super::FaceWindingHonesty]| -> Vec<_> {
        let mut k: Vec<_> = vs.iter().map(|v| (v.face, v.honest, v.reason)).collect();
        k.sort_by_key(|x| x.0);
        k
    };
    assert_eq!(
        key(&among),
        key(&whole),
        "{label}: whole-region audit differs"
    );
}

#[test]
fn region_audit_matches_the_whole_body_audit_on_primitives() {
    for (label, b) in bodies() {
        assert_region_matches_whole(label, &b);
    }
}

#[test]
fn region_audit_follows_tangency_witness_chains() {
    let b = chained_torus(0, None);
    let whole = b.audit_winding_honesty();
    let chained = whole
        .iter()
        .filter(|v| v.reason == "tangency-witness")
        .count();
    assert!(
        chained + 1 >= whole.len(),
        "the fixture must resolve through witness chains ({chained}/{})",
        whole.len()
    );
    assert_region_matches_whole("chained torus", &b);
    // A dishonest link inside the chain: its own and its dependants'
    // verdicts must still agree with the whole-body audit.
    for flip in [1, 5, 10] {
        assert_region_matches_whole(
            "chained torus, one flag flipped",
            &chained_torus(0, Some(flip)),
        );
    }
}

#[test]
fn region_reconcile_flips_only_the_region() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces: Vec<FaceHandle> = cube.faces.keys().collect();
    for &f in &faces[..2] {
        let s = cube.faces[f].surface.as_mut().unwrap();
        s.orientation = !s.orientation;
    }
    assert_eq!(cube.reconcile_orientation_flags_among(&faces[1..]), 1);
    let dishonest = cube.dishonest_winding_faces();
    assert_eq!(
        dishonest,
        vec![faces[0]],
        "the face outside the region keeps its flag"
    );
    assert_eq!(cube.reconcile_orientation_flags_among(&faces[..1]), 1);
    assert!(cube.dishonest_winding_faces().is_empty());
}
