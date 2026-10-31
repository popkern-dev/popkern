//! The face-side rule (plan §4.3a gate): a group growing across a
//! boundary edge of a face does not take points on both sides of it.

use popkern_core::geometry::Point3;
use popkern_core::topology::FaceHandle;

use super::candidates::{CandidateVertex, Incidence, SourceKind};
use super::resolve::{Incompatibility, ResolutionContext, resolve};
use super::tests::{TOL, cube, options, record_end};

/// A point of `face` (operand 0) at `at`.
fn on_face(face: FaceHandle, at: Point3, key: u64) -> CandidateVertex {
    CandidateVertex {
        position: at,
        radius: TOL,
        zone: TOL,
        source: SourceKind::EdgeFaceHit,
        incidence: vec![Incidence::Face { operand: 0, face }],
        key: (SourceKind::EdgeFaceHit, 1, key, 0, 0),
    }
}

/// Two points of the cube's top face 1.5 tol either side of its x = 1
/// boundary edge, linked through a record end on the edge: the chain
/// would grow into one vertex 1.5 tol wide, but the two face points
/// lie on opposite sides of the edge, so they stay apart.
#[test]
fn a_group_does_not_span_a_boundary_edge_of_a_face_both_sides_lie_on() {
    let c = cube();
    let top = c
        .faces
        .keys()
        .find(|&f| {
            let ps = c.face_outer_loop_positions(f);
            !ps.is_empty() && ps.iter().all(|p| (p.z - 1.0).abs() < 1e-9)
        })
        .expect("top face");
    let cands = vec![
        on_face(top, Point3::new(1.0 - 1.5 * TOL, 0.5, 1.0), 0),
        record_end(0, 0.0, false, Point3::new(1.0, 0.5, 1.0)),
        on_face(top, Point3::new(1.0 + 1.5 * TOL, 0.5, 1.0), 1),
    ];
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &[],
    };
    let res = resolve(&cands, &ctx, options(10.0 * TOL)).expect("resolves");
    assert_ne!(
        res.group_of[0], res.group_of[2],
        "both sides in one vertex: {res:?}"
    );
    assert!(
        res.separations
            .iter()
            .any(|&(_, _, why)| why == Incompatibility::OppositeFaceSides),
        "{:?}",
        res.separations
    );

    // Both on the same side, 1.5 and 0.5 tol in from the edge: one vertex.
    let same_side = vec![
        on_face(top, Point3::new(1.0 - 1.5 * TOL, 0.5, 1.0), 0),
        record_end(0, 0.0, false, Point3::new(1.0, 0.5, 1.0)),
        on_face(top, Point3::new(1.0 - 0.5 * TOL, 0.5, 1.0), 1),
    ];
    let res = resolve(&same_side, &ctx, options(10.0 * TOL)).expect("resolves");
    assert_eq!(res.groups.len(), 1, "{res:?}");
}
