//! Blends refuse a body they do not accept by its kind (PK-1101 S5),
//! before any geometry is read: edge fillets and chamfers, their edits
//! and vertex blends take solids, corner blends sheets or solids.

use popkern_core::error::CoreError;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

use crate::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// A cube beside a sheet: lumps of two kinds, a general body.
fn general_body() -> BRep {
    let mut body = BRep::unit_cube(tol()).expect("cube");
    let mut sheet = BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(3.0, 0.0, 0.0)))
        .expect("cube");
    let face = sheet.faces.keys().next().expect("a face");
    sheet.retire_faces_to_boundary(&[face]).expect("retires");
    sheet.settle_open_shells();
    body.append_body(&sheet);
    assert_eq!(body.body_kind(), BodyKind::General);
    body
}

fn assert_refused<T: std::fmt::Debug>(result: Result<T, BlendError>, operation: &str) {
    match result {
        Err(BlendError::Topology(CoreError::UnsupportedBodyKind {
            operation: op,
            found: BodyKind::General,
            ..
        })) => assert_eq!(op, operation),
        other => panic!("{operation}: expected the kind gate, got {other:?}"),
    }
}

#[test]
fn blends_refuse_a_general_body() {
    let body = general_body();
    let edge = body.edges.keys().next().expect("an edge");
    let vertex = body.vertices.keys().next().expect("a vertex");
    let face = body.faces.keys().next().expect("a face");
    let fillet = [FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let chamfer = [ChamferSpec {
        edge,
        kind: ChamferKind::Equal(0.1),
    }];
    assert_refused(
        fillet_edges(&body, &fillet, &FilletOptions::default(), tol(), None),
        "fillet",
    );
    assert_refused(
        chamfer_edges(&body, &chamfer, &ChamferOptions::default(), tol(), None),
        "chamfer",
    );
    assert_refused(fillet_edge(&body, edge, 0.1, tol(), None), "fillet");
    assert_refused(fillet_corner_2d(&body, vertex, 0.1), "corner fillet");
    assert_refused(chamfer_corner_2d(&body, vertex, 0.1, 0.1), "corner chamfer");
    assert_refused(
        modify_fillet_radius(&body, face, 0.2, tol()),
        "modify fillet",
    );
    assert_refused(
        modify_chamfer_distance(&body, face, ChamferKind::Equal(0.2), tol()),
        "modify chamfer",
    );
    assert_refused(unfillet_face(&body, face), "remove fillet");
    assert_refused(
        vertex_blend_3d(&body, vertex, 0.1, &FilletOptions::default(), tol(), None),
        "vertex blend",
    );
}
