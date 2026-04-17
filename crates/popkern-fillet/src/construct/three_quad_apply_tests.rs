//! Integration tests for the three-quad apply pipeline:
//! * R8 C6 regression guard — mixed-radius cube corner lands as a
//!   manifold hex (6 HE, 6 distinct vertices, NURBS surface).
//! * R8 C7 dispatch — when `use_rational_gregory` is on, the
//!   mixed-radius corner's blend face is poked into three
//!   rectangular-Gregory sub-patches.

use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::{CornerMethod, fillet_edges};

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn r8_c6_mixed_radius_cube_corner_lands_manifold_hex() {
    // Fixture: mixed-radius cube corner with 3 fillets. The R8 C6
    // fix (per-face trim anchors + kev-for-bridge) must produce a
    // 6-half-edge blend face with 6 distinct boundary vertices.
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    assert_eq!(incident_edges.len(), 3);

    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let tol = Tolerance::default();
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol, None).unwrap();

    // Locate the corner patch for our target vertex.
    let target_info = result
        .corner_patches
        .iter()
        .find(|p| p.vertex == vertex)
        .expect("target vertex should have a corner patch");
    let hex_face = target_info.face;
    assert_eq!(
        target_info.method,
        CornerMethod::Gregory,
        "mixed-radius corner should route to Gregory (polar-NURBS) in R8 C6",
    );

    // Regression guard: the blend face must survive the post-insert
    // pipeline. Before R8 C6 this face (a pinched pentagon) was
    // getting absorbed by `sliver_absorb::find_redundant_edge_target`
    // because `surface_fix` reassigned its NURBS to a Cylinder
    // matching an adjacent fillet strip.
    let brep = &result.brep;
    assert!(
        brep.faces.contains_key(hex_face),
        "R8 C6 regression: corner blend face must not be absorbed by sliver_absorb"
    );

    // Walk the face's outer loop, collect vertex handles.
    let loop_handle = brep.faces[hex_face].outer_loop;
    let start_he = brep.loops[loop_handle].half_edge;
    let mut count = 0usize;
    let mut cursor = start_he;
    let mut verts = Vec::new();
    loop {
        count += 1;
        if count > 12 {
            panic!("blend-face loop doesn't close (>12 half-edges)");
        }
        verts.push(brep.half_edges[cursor].origin);
        cursor = brep.half_edges[cursor].next;
        if cursor == start_he {
            break;
        }
    }

    assert_eq!(
        count, 6,
        "R8 C6: mixed-radius corner blend face must be a hexagon; \
         got {count} half-edges, vertices {verts:?}"
    );

    // Distinctness: no vertex appears twice (no pinch). Before
    // R8 C6 the pinched pentagon had one vertex visited twice.
    let mut unique = verts.clone();
    unique.sort();
    unique.dedup();
    assert_eq!(
        unique.len(),
        6,
        "R8 C6: hex should have 6 distinct boundary vertices; \
         got {} unique from {verts:?}",
        unique.len()
    );
}

/// R8 C7 dispatch test: with `use_rational_gregory` enabled, the
/// mixed-radius cube corner's hex blend face is poked into three
/// rectangular Gregory sub-patches — +1 central vertex, +3 spoke
/// edges, +2 sub-faces (Euler invariant `ΔV - ΔE + ΔF = 0`).
#[test]
fn r8_c7_dispatch_produces_three_rectangular_gregory_subpatches() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    assert_eq!(incident_edges.len(), 3);

    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();

    // Baseline (dispatch off): single corner patch, NURBS surface.
    let baseline = fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .unwrap();
    let baseline_corners: Vec<_> = baseline
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .collect();
    assert_eq!(
        baseline_corners.len(),
        1,
        "baseline: a single corner patch per vertex"
    );

    // Rational-Gregory dispatch on.
    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None).unwrap();

    let gregory_corners: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .collect();
    assert_eq!(
        gregory_corners.len(),
        3,
        "R8 C7: rational-Gregory dispatch should yield 3 corner patches per \
         mixed-radius trihedral vertex (one per quad sub-patch)"
    );

    for patch in &gregory_corners {
        assert_eq!(
            patch.method,
            CornerMethod::Gregory,
            "R8 C7: sub-patches report Gregory as their method"
        );
        let sof = result.brep.faces[patch.face]
            .surface
            .as_ref()
            .expect("R8 C7: sub-face must have a surface");
        match &sof.geometry {
            SurfaceGeometry::Gregory(GregoryPatch::Rectangular(_)) => {}
            other => panic!(
                "R8 C7: sub-face surface must be a rectangular Gregory; \
                 got {other:?}"
            ),
        }
    }

    // Euler invariant: +1V, +3E, +2F compared with the baseline.
    assert_eq!(
        result.brep.vertices.len(),
        baseline.brep.vertices.len() + 1,
        "R8 C7: central vertex added"
    );
    assert_eq!(
        result.brep.edges.len(),
        baseline.brep.edges.len() + 3,
        "R8 C7: three spoke edges added"
    );
    assert_eq!(
        result.brep.faces.len(),
        baseline.brep.faces.len() + 2,
        "R8 C7: hex splits into three sub-faces (+2 faces)"
    );
}
