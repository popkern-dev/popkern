//! Integration test for the three-quad apply helper. Runs the
//! full `fillet_edges` pipeline on a mixed-radius cube corner,
//! locates the corner-blend face in the result, and verifies the
//! R8 C6 hex topology invariants.
//!
//! The three-quad *apply* wiring (poke + Gregory assignment) is
//! R8 C7 follow-up — it needs the blend boundary cached at insert
//! time so it survives the R8 C6 `kev` that removes V.

use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::{CornerMethod, fillet_edges};

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
