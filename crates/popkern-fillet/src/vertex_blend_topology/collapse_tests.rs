//! Tests for `collapse_apex_lens`.
//!
//! Two test categories:
//! 1. **Precondition rejection** — cube_uniform's lens (cross-
//!    cylinder junction) should be rejected with a structured
//!    `CollapseError::CrossFaceJunction`.
//! 2. **Happy path** — a synthetically-built BRep with a
//!    same-twin-face lens should be collapsed cleanly, leaving
//!    a valid BRep with the expected vertex / edge / face deltas.

use super::collapse::{CollapseError, collapse_apex_lens};
use super::execute::execute_deletion_plan;
use super::plan::plan_apex_region_deletion;
use crate::FilletOptions;
use crate::types::{FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

fn corner_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    use std::collections::BTreeSet;
    brep.half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn cube_apex(brep: &BRep) -> VertexHandle {
    brep.vertices
        .keys()
        .find(|&v| brep.half_edges.values().filter(|he| he.origin == v).count() == 3)
        .expect("cube has at least one valence-3 vertex")
}

fn fillet_cube_uniform_with_lenses() -> (BRep, Vec<FaceHandle>, Vec<HalfEdgeHandle>) {
    // Run the fillet pipeline AND Phase 6c-C.2b's mef step so
    // we have lens faces to feed to the collapse op.
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let apex_pre = cube_apex(&cube);
    let edges = corner_edges(&cube, apex_pre);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.10),
        })
        .collect();
    let result = crate::fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("fillet uniform cube");
    let cap_face = result.corner_patches[0].face;
    let plan = plan_apex_region_deletion(&result.brep, cap_face);
    let mut brep = result.brep;
    let exec = execute_deletion_plan(&mut brep, &plan).expect("mef step ok");
    // For each lens face, find its 2 generator he's in the
    // lens loop. The third he is the new arc inserted by mef.
    let mut gen_pairs: Vec<HalfEdgeHandle> = Vec::new();
    for &lens in &exec.lens_faces {
        let lens_loop = brep.faces[lens].outer_loop;
        let start = brep.loops[lens_loop].half_edge;
        let mut cursor = start;
        for _ in 0..3 {
            let he_edge = brep.half_edges[cursor].edge;
            let he_next = brep.half_edges[cursor].next;
            if !exec.inserted_arcs.contains(&he_edge) {
                gen_pairs.push(cursor);
            }
            cursor = he_next;
            if cursor == start {
                break;
            }
        }
    }
    (brep, exec.lens_faces, gen_pairs)
}

#[test]
#[ignore = "route B corner reconstruction (contact collapse + exact boundary arcs) removed the legacy anchor-hexagon cap anatomy this machinery was built to clean up; the apex-region deletion plan finds no redundant lens structures on the reconstructed corner. Retire or repurpose with PK-101."]
fn cube_uniform_collapse_all_three_lenses() {
    // Run mef + collapse for all 3 lenses. Expected post-state:
    // cap face has 3 outer-loop edges (the cap-tangent arcs);
    // 3 footprint vertices removed; 6 generator edges removed;
    // 3 lens faces removed.
    let (mut brep, lens_faces, gens) = fillet_cube_uniform_with_lenses();
    let v_before = brep.vertices.len();
    let e_before = brep.edges.len();
    let f_before = brep.faces.len();
    for i in 0..3 {
        let lens = lens_faces[i];
        let gen1 = gens[2 * i];
        let gen2 = gens[2 * i + 1];
        let _ =
            collapse_apex_lens(&mut brep, lens, gen1, gen2).expect("each lens collapses cleanly");
    }
    assert_eq!(
        brep.vertices.len(),
        v_before - 3,
        "3 footprint vertices removed"
    );
    assert_eq!(brep.edges.len(), e_before - 6, "6 generator edges removed");
    assert_eq!(brep.faces.len(), f_before - 3, "3 lens faces removed");
    brep.validate()
        .expect("BRep validates after all 3 collapses");
}

#[test]
#[ignore = "route B corner reconstruction (contact collapse + exact boundary arcs) removed the legacy anchor-hexagon cap anatomy this machinery was built to clean up; the apex-region deletion plan finds no redundant lens structures on the reconstructed corner. Retire or repurpose with PK-101."]
fn cube_uniform_lens_collapse_diagnostic() {
    // Diagnostic — figure out whether cube_uniform's lens IS
    // actually a same-twin-face case (in which case we can
    // collapse it directly) or a cross-cylinder junction
    // (which the simple primitive must reject).
    use popkern_core::geometry::surface::SurfaceGeometry;
    let (mut brep, lens_faces, gens) = fillet_cube_uniform_with_lenses();
    assert_eq!(lens_faces.len(), 3);
    assert_eq!(gens.len(), 6);

    let lens0 = lens_faces[0];
    let gen1 = gens[0];
    let gen2 = gens[1];

    let gen1_twin = brep.half_edges[gen1].twin;
    let gen2_twin = brep.half_edges[gen2].twin;
    let f1 = brep.loops[brep.half_edges[gen1_twin].loop_].face;
    let f2 = brep.loops[brep.half_edges[gen2_twin].loop_].face;
    let surf_kind = |f: popkern_core::topology::FaceHandle| -> &'static str {
        brep.faces[f]
            .surface
            .as_ref()
            .map(|s| match &s.geometry {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cyl",
                SurfaceGeometry::Sphere { .. } => "Sph",
                _ => "Other",
            })
            .unwrap_or("?")
    };
    eprintln!(
        "lens0={lens0:?} gen1.twin.face={f1:?} ({}) gen2.twin.face={f2:?} ({})",
        surf_kind(f1),
        surf_kind(f2),
    );
    let snapshot = brep.clone();
    let result = collapse_apex_lens(&mut brep, lens0, gen1, gen2);
    eprintln!("collapse result: {result:?}");
    if result.is_ok() {
        // Verify the post-state is still valid and report the
        // BRep deltas.
        eprintln!(
            "delta: vertices {} -> {}, edges {} -> {}, faces {} -> {}",
            snapshot.vertices.len(),
            brep.vertices.len(),
            snapshot.edges.len(),
            brep.edges.len(),
            snapshot.faces.len(),
            brep.faces.len(),
        );
    }
}

#[test]
#[ignore = "route B corner reconstruction (contact collapse + exact boundary arcs) removed the legacy anchor-hexagon cap anatomy this machinery was built to clean up; the apex-region deletion plan finds no redundant lens structures on the reconstructed corner. Retire or repurpose with PK-101."]
fn rejects_invalid_lens_face_handle() {
    let (mut brep, _, _) = fillet_cube_uniform_with_lenses();
    // Use an arbitrary half-edge handle as a fake lens face —
    // the function should reject because the face handle is
    // invalid.
    let any_he = brep.half_edges.keys().next().unwrap();
    // Construct an obviously-invalid FaceHandle by reusing
    // a slot that's been removed.  Easier: use a real face
    // but pass random half-edge handles for gen1/gen2.
    let any_face = brep.faces.keys().next().unwrap();
    let err = collapse_apex_lens(&mut brep, any_face, any_he, any_he)
        .expect_err("invalid input rejected");
    assert!(matches!(
        err,
        CollapseError::LensNotTriangle { .. }
            | CollapseError::GenNotInLens
            | CollapseError::DuplicateGen
            | CollapseError::CrossFaceJunction
            | CollapseError::TwinNotConsecutive
            | CollapseError::GensDontMeetAtSharedVertex
            | CollapseError::ThirdEdgeFolds
    ));
}
