//! PK-194 — control cage to smooth solid, through the Part API.
//!
//! These gates prove the polyhedral-net spline family is reachable as
//! a product capability rather than a core-internal experiment: a
//! caller hands over a control cage and gets back a validated `Part`
//! that the rest of the kernel treats like any other body.
//!
//! The cages here are deliberately *not* the TOMS example fixtures —
//! the point is that the path works on geometry written for this test,
//! not only on the inputs the family was developed against.

use popkern_core::geometry::Point3;
use popkern_core::geometry::polyhedral_spline::ControlNet;
use popkern_core::tolerance::Tolerance;
use popkern_part::{SmoothMeshOptions, smooth_solid_from_control_net, smooth_solid_from_obj};

/// A closed cube cage: 8 vertices, 6 quads, every vertex valence 3.
///
/// The regular-quad interior of this cage is the simplest closed
/// control net the family can consume, and every vertex is an
/// extraordinary point (valence 3 ≠ 4), so it exercises the
/// `ExtraordinaryStar` construction at every corner rather than only
/// the regular case.
fn cube_cage() -> &'static str {
    "\
# unit cube control cage
v -1 -1 -1
v  1 -1 -1
v  1  1 -1
v -1  1 -1
v -1 -1  1
v  1 -1  1
v  1  1  1
v -1  1  1
f 1 4 3 2
f 5 6 7 8
f 1 2 6 5
f 2 3 7 6
f 3 4 8 7
f 4 1 5 8
"
}

/// A cage with a valence-5 vertex: the cube with one corner truncated
/// into a pentagon, so the net contains an n-gon face and a
/// higher-valence star.
fn truncated_corner_cage() -> ControlNet {
    // Cube corners with corner (1,1,1) cut off by a triangle through
    // the three adjacent edge midpoints.
    let vertices = vec![
        Point3::new(-1.0, -1.0, -1.0),
        Point3::new(1.0, -1.0, -1.0),
        Point3::new(1.0, 1.0, -1.0),
        Point3::new(-1.0, 1.0, -1.0),
        Point3::new(-1.0, -1.0, 1.0),
        Point3::new(1.0, -1.0, 1.0),
        Point3::new(-1.0, 1.0, 1.0),
        // truncation triangle
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 1.0),
    ];
    let faces = vec![
        vec![0, 3, 2, 1],
        vec![4, 5, 7, 9, 6],
        vec![0, 1, 5, 4],
        vec![1, 2, 8, 7, 5],
        vec![2, 3, 6, 9, 8],
        vec![3, 0, 4, 6],
        vec![7, 8, 9],
    ];
    ControlNet::new(vertices, faces).expect("truncated-corner cage must be a valid control net")
}

#[test]
fn cube_cage_becomes_a_smooth_part() {
    let (part, report) = smooth_solid_from_obj(
        cube_cage(),
        Tolerance::default(),
        SmoothMeshOptions::default(),
    )
    .expect("a closed cube cage must produce a smooth solid");

    assert_eq!(report.source_vertices, 8, "cage vertices consumed");
    assert_eq!(report.source_faces, 6, "cage faces consumed");
    assert!(
        report.accepted_components > 0,
        "at least one spline component must be realised"
    );
    assert_eq!(
        report.remaining_boundary_edges, 0,
        "a closed cage must weld into a closed body, {} edge(s) left over",
        report.remaining_boundary_edges
    );
    assert!(
        part.brep().validate().is_ok(),
        "the smooth body must pass kernel validation"
    );
    assert!(
        part.brep().faces.len() > 6,
        "the spline body should refine the 6-face cage, got {} faces",
        part.brep().faces.len()
    );
}

#[test]
fn smooth_body_is_ordinary_geometry_downstream() {
    // The point of returning a Part: the body is not a special case.
    // Tessellation is the cheapest proof that downstream consumers
    // accept it.
    let (part, _) = smooth_solid_from_obj(
        cube_cage(),
        Tolerance::default(),
        SmoothMeshOptions::default(),
    )
    .expect("cube cage must build");
    let mesh = popkern_mesh::tessellate(part.brep(), &popkern_mesh::TessellationOptions::default())
        .expect("the smooth body must tessellate");
    assert!(
        mesh.triangles.len() > 12,
        "a smooth body should tessellate to more than the cage's own facets, got {}",
        mesh.triangles.len()
    );
}

#[test]
fn n_gon_and_high_valence_cage_is_served_or_declines_with_reasons() {
    let net = truncated_corner_cage();
    let result =
        smooth_solid_from_control_net(&net, Tolerance::default(), SmoothMeshOptions::default());
    match result {
        Ok((part, report)) => {
            assert!(
                part.brep().validate().is_ok(),
                "an accepted n-gon cage must still validate"
            );
            assert_eq!(
                report.remaining_boundary_edges, 0,
                "accepted body is closed"
            );
        }
        Err(e) => {
            // A decline is acceptable for this harder cage, but it must
            // say which stage failed and why — a silent or generic
            // failure would leave the caller with nothing to act on.
            let message = e.to_string();
            assert!(
                message.contains("smooth mesh failed at"),
                "decline must name the pipeline stage, got: {message}"
            );
            assert!(
                message.len() > 40,
                "decline must carry a real reason, got: {message}"
            );
        }
    }
}

#[test]
fn open_cage_declines_unless_sheet_bodies_are_allowed() {
    // Five faces of the cube: deliberately open.
    let open = "\
v -1 -1 -1
v  1 -1 -1
v  1  1 -1
v -1  1 -1
v -1 -1  1
v  1 -1  1
v  1  1  1
v -1  1  1
f 1 4 3 2
f 1 2 6 5
f 2 3 7 6
f 3 4 8 7
f 4 1 5 8
";
    let closed_required = smooth_solid_from_obj(
        open,
        Tolerance::default(),
        SmoothMeshOptions {
            require_closed: true,
            ..SmoothMeshOptions::default()
        },
    );
    match closed_required {
        Err(e) => {
            let message = e.to_string();
            assert!(
                message.contains("smooth mesh failed at"),
                "open cage must decline with a staged reason, got: {message}"
            );
        }
        Ok((_, report)) => panic!(
            "an open cage must not silently pass as a closed solid (report left {} boundary edges)",
            report.remaining_boundary_edges
        ),
    }
}

#[test]
fn malformed_obj_declines_at_the_parse_stage() {
    let bad = "v 0 0\nf 1 2 3\n";
    let err = smooth_solid_from_obj(bad, Tolerance::default(), SmoothMeshOptions::default())
        .err()
        .expect("a malformed cage must not produce a body");
    let message = err.to_string();
    assert!(
        message.contains("obj parse"),
        "parse failure must be attributed to the parse stage, got: {message}"
    );
}
