//! Tests that shell operations produce pcurves on all face half-edges.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn shelled_box_has_pcurves_on_all_half_edges() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[first_face],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // Check half-edges reachable from face loops (assembly may leave
    // orphan half-edges from removed-face edges in the slotmap).
    for (_fh, face) in result.brep.faces.iter() {
        let mut loops_to_check = vec![face.outer_loop];
        loops_to_check.extend_from_slice(&face.inner_loops);
        for lh in loops_to_check {
            if let Some(loop_) = result.brep.loops.get(lh) {
                let start = loop_.half_edge;
                let mut cur = start;
                loop {
                    let he = &result.brep.half_edges[cur];
                    assert!(
                        he.pcurve.is_some(),
                        "shell half-edge {:?} missing pcurve",
                        cur
                    );
                    cur = he.next;
                    if cur == start {
                        break;
                    }
                }
            }
        }
    }
}

#[test]
fn shell_no_removal_has_pcurves() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    for (_fh, face) in result.brep.faces.iter() {
        let mut loops_to_check = vec![face.outer_loop];
        loops_to_check.extend_from_slice(&face.inner_loops);
        for lh in loops_to_check {
            if let Some(loop_) = result.brep.loops.get(lh) {
                let start = loop_.half_edge;
                let mut cur = start;
                loop {
                    let he = &result.brep.half_edges[cur];
                    assert!(
                        he.pcurve.is_some(),
                        "shell (no removal) half-edge {:?} missing pcurve",
                        cur
                    );
                    cur = he.next;
                    if cur == start {
                        break;
                    }
                }
            }
        }
    }
}
