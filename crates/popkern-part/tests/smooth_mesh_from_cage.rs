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
    if !toms_tables_or_skip() {
        return;
    }
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
    // PK-198 no-regression pin: the all-quad-neighbourhood cage keeps
    // its direct serve — no cage normalization, machine-zero welds.
    assert_eq!(
        report.refinement_steps, 0,
        "the cube cage is family-native and must not be refined"
    );
    assert!(report.attempts.is_empty(), "no failed attempts");
    assert!(
        report.max_weld_gap < 1e-12,
        "cube welds are machine-zero, got {:.3e}",
        report.max_weld_gap
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
    if !toms_tables_or_skip() {
        return;
    }
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
fn n_gon_and_high_valence_cage_serves_through_cage_normalization() {
    if !toms_tables_or_skip() {
        return;
    }
    // PK-198: the direct enumeration covers exactly one vertex of this
    // cage (the only all-quad neighbourhood), so the family cannot
    // close it as given. One Catmull–Clark step normalizes it into an
    // all-quad cage every vertex of which is served, and the body
    // closes and validates.
    let net = truncated_corner_cage();
    let (part, report) =
        smooth_solid_from_control_net(&net, Tolerance::default(), SmoothMeshOptions::default())
            .expect("the truncated-corner cage must serve through normalization");
    assert_eq!(report.refinement_steps, 1, "one Catmull–Clark step");
    assert_eq!(report.attempts.len(), 1, "the direct attempt is recorded");
    assert!(
        report.attempts[0].contains("unserved"),
        "the recorded direct failure explains the coverage gap, got: {}",
        report.attempts[0]
    );
    assert_eq!(report.remaining_boundary_edges, 0, "closed body");
    assert_eq!(report.declined_pairs, 0, "every boundary pair welds");
    assert!(report.rejections.is_empty(), "no spline rejections");
    assert!(part.brep().validate().is_ok(), "kernel validation green");
    assert!(
        part.brep().validate_all().is_empty(),
        "full validation (incl. self-intersection) green"
    );
}

#[test]
fn concave_hexagonal_cage_serves_through_cage_normalization() {
    if !toms_tables_or_skip() {
        return;
    }
    // PK-198 l_block class: hexagon end caps put a valence-6 face
    // point into the refined cage, whose star tables tile at
    // half-quarter granularity — the complex-level tiling unification
    // subdivides the coarser components so the weld pairs 1:1.
    let l_block = "\
v 0 0 0
v 2 0 0
v 2 1 0
v 1 1 0
v 1 2 0
v 0 2 0
v 0 0 1
v 2 0 1
v 2 1 1
v 1 1 1
v 1 2 1
v 0 2 1
f 6 5 4 3 2 1
f 7 8 9 10 11 12
f 1 2 8 7
f 2 3 9 8
f 3 4 10 9
f 4 5 11 10
f 5 6 12 11
f 6 1 7 12
";
    let (part, report) =
        smooth_solid_from_obj(l_block, Tolerance::default(), SmoothMeshOptions::default())
            .expect("the concave l_block cage must serve through normalization");
    assert_eq!(report.refinement_steps, 1, "one Catmull–Clark step");
    assert_eq!(report.remaining_boundary_edges, 0, "closed body");
    assert_eq!(report.declined_pairs, 0, "every boundary pair welds");
    assert!(part.brep().validate().is_ok(), "kernel validation green");
    assert!(
        part.brep().validate_all().is_empty(),
        "full validation (incl. self-intersection) green"
    );
}

#[test]
fn unsupported_valence_declines_naming_the_vertex() {
    // A 9-sided pyramid: the apex has valence 9, above the family's
    // 3..=8 window. Refinement preserves valence, so the decline is
    // structural and must name the offending vertex.
    let mut vertices = vec![Point3::new(0.0, 0.0, 1.0)];
    for k in 0..9 {
        let angle = 2.0 * std::f64::consts::PI * k as f64 / 9.0;
        vertices.push(Point3::new(angle.cos(), angle.sin(), 0.0));
    }
    let mut faces: Vec<Vec<usize>> = (0..9).map(|k| vec![0, 1 + k, 1 + (k + 1) % 9]).collect();
    faces.push((1..=9).rev().collect());
    let net = ControlNet::new(vertices, faces).expect("pyramid cage");
    let err =
        smooth_solid_from_control_net(&net, Tolerance::default(), SmoothMeshOptions::default())
            .expect_err("valence 9 is outside the family's window");
    let message = err.to_string();
    assert!(
        message.contains("vertex 0") && message.contains("valence 9"),
        "decline must name the apex and its valence, got: {message}"
    );
}

#[test]
fn non_manifold_cage_declines_naming_the_edge() {
    // Three faces share the edge 0-1: the halfedge construction
    // reports the duplicated directed edge by name before any
    // geometry is built.
    let net = ControlNet::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(0.5, -1.0, 0.0),
            Point3::new(0.5, 0.0, 1.0),
        ],
        vec![vec![0, 1, 2], vec![1, 0, 3], vec![0, 1, 4]],
    )
    .expect("indices are valid even when non-manifold");
    let err =
        smooth_solid_from_control_net(&net, Tolerance::default(), SmoothMeshOptions::default())
            .expect_err("a non-manifold cage must decline");
    let message = err.to_string();
    assert!(
        message.contains("duplicate directed edge 0->1"),
        "decline must name the non-manifold edge, got: {message}"
    );
}

#[test]
fn self_intersecting_body_declines_naming_the_faces() {
    if !toms_tables_or_skip() {
        return;
    }
    // A cube cage with one corner pulled just past the opposite
    // corner: the direct smooth body grazes itself, and the
    // validate_all screen declines naming the face pair. (Deeper
    // pulls contract into embedded — measured crossing-free — bodies;
    // this marginal graze is the reaching case.)
    let mut vertices = vec![
        Point3::new(-1.0, -1.0, -1.0),
        Point3::new(1.0, -1.0, -1.0),
        Point3::new(1.0, 1.0, -1.0),
        Point3::new(-1.0, 1.0, -1.0),
        Point3::new(-1.0, -1.0, 1.0),
        Point3::new(1.0, -1.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(-1.0, 1.0, 1.0),
    ];
    vertices[6] = Point3::new(-2.0, -2.0, -2.0);
    let faces = vec![
        vec![0, 3, 2, 1],
        vec![4, 5, 6, 7],
        vec![0, 1, 5, 4],
        vec![1, 2, 6, 5],
        vec![2, 3, 7, 6],
        vec![3, 0, 4, 7],
    ];
    let net = ControlNet::new(vertices, faces).expect("pulled-corner cage");
    let err = smooth_solid_from_control_net(
        &net,
        Tolerance::default(),
        SmoothMeshOptions {
            max_refinement_steps: 0,
            ..SmoothMeshOptions::default()
        },
    )
    .expect_err("the grazing body must decline the self-intersection screen");
    let message = err.to_string();
    assert!(
        message.contains("self-intersection") && message.contains("face"),
        "decline must name the screen and the faces, got: {message}"
    );

    // With normalization allowed, the refined realization of the same
    // cage is tighter and measured crossing-free — the retry serves a
    // genuinely clean body rather than masking the decline.
    let (part, report) =
        smooth_solid_from_control_net(&net, Tolerance::default(), SmoothMeshOptions::default())
            .expect("the refined realization serves clean");
    assert_eq!(report.refinement_steps, 1);
    assert!(part.brep().validate_all().is_empty());
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
        .expect_err("a malformed cage must not produce a body");
    let message = err.to_string();
    assert!(
        message.contains("obj parse"),
        "parse failure must be attributed to the parse stage, got: {message}"
    );
}

/// The PnS mask tables are optional (ACM ToMS package, untracked —
/// PK-314/PK-315). Tests that build smooth-mesh bodies need them at
/// runtime; on checkouts without the package they skip with a visible
/// notice instead of failing.
fn toms_tables_or_skip() -> bool {
    if popkern_core::geometry::polyhedral_spline::toms_tables_available() {
        return true;
    }
    eprintln!("[toms] ACM ToMS 1032 package absent — smooth-mesh test skipped (PK-314/PK-315)");
    false
}
