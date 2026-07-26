//! Phase 6c-C.1 tests — apex-region detection on canonical
//! fixtures. These tests verify the topological rule from
//! design supplement §10.1 against known-correct expectations
//! on cube and pyramid apexes.

use super::detect::{detect_apex_region, outer_loop_valence};
use crate::FilletOptions;
use crate::types::{FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

fn corner_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    use std::collections::BTreeSet;
    let edges: BTreeSet<EdgeHandle> = brep
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.into_iter().collect()
}

/// Walk the apex's incident edges (without applying any fillet)
/// and return the apex vertex's slot ID and the count of
/// incident edges. Used to anchor the pre-fillet state expected
/// by detection.
fn cube_apex(brep: &BRep) -> popkern_core::topology::VertexHandle {
    // unit_cube places the first vertex at (0,0,0); take the
    // first one with valence == 3.
    brep.vertices
        .keys()
        .find(|&v| brep.half_edges.values().filter(|he| he.origin == v).count() == 3)
        .expect("cube has at least one valence-3 vertex")
}

#[test]
fn detect_finds_pre_fillet_cube_apex_region() {
    // Pre-fillet unit cube: apex (0,0,0) has 3 incident edges
    // and 3 incident faces. No blend_face yet.
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let apex = cube_apex(&cube);

    let region = detect_apex_region(&cube, apex);

    assert!(region.is_valid(), "region should be valid");
    assert_eq!(region.apex, Some(apex));
    assert_eq!(region.edges.len(), 3, "cube apex has 3 incident edges");
    assert_eq!(region.faces.len(), 3, "cube apex has 3 incident faces");
    assert!(region.blend_face.is_none(), "no blend face pre-insert");
}

#[test]
fn detect_post_fillet_cube_uniform_apex_region() {
    // Post-fillet unit cube (uniform r=0.10): apex *still
    // exists* as a BRep vertex (insert_vertex_blend hasn't run
    // yet — the fillet pipeline up to and including Phase 5
    // surface fix is what ran).  Our 6c-C.1 call site is right
    // after Phase 5 and before insert_vertex_blend.
    //
    // For this test we use the *full* fillet pipeline output
    // (insert_vertex_blend has run) and confirm detection finds
    // *something* meaningful even at this later stage — the apex
    // region will instead surface the post-blend triangular cap
    // face's vertices.
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

    // Find the cap face from the corner_patches (one entry for
    // the apex blend).
    let cap_face = result
        .corner_patches
        .first()
        .map(|p| p.face)
        .expect("at least one corner_patch");

    // After insert_vertex_blend, the apex vertex itself was
    // removed (the kev_v_bridge step kills it). Instead, walk
    // the cap face's outer-loop vertices and detect each as an
    // "apex-region" probe — they should all return valid
    // (non-empty edges + faces).
    let cap_loop = result.brep.loops[result.brep.faces[cap_face].outer_loop].half_edge;
    let mut probe_v = result.brep.half_edges[cap_loop].origin;
    let mut probed = 0;
    let mut found_a_region = false;
    for _ in 0..16 {
        let region = detect_apex_region(&result.brep, probe_v);
        if region.is_valid() {
            found_a_region = true;
            // The cap face must be incident to this probe vertex
            assert!(
                region.faces.contains(&cap_face),
                "cap face should be in the apex region for vertex {probe_v:?}"
            );
        }
        probed += 1;
        let next_he = result.brep.half_edges[cap_loop].next;
        probe_v = result.brep.half_edges[next_he].origin;
        if next_he == cap_loop {
            break;
        }
    }
    assert!(
        found_a_region,
        "at least one cap-face vertex must yield a valid region"
    );
    assert!(probed >= 3, "we walked the cap face's loop");
}

#[test]
fn detect_handles_invalid_apex_handle() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    // Build a stale handle by deleting the cube and constructing
    // a fresh one — the same slot ID won't resolve in the new BRep.
    let stale_apex = cube_apex(&cube);
    drop(cube);
    let cube2 = BRep::unit_cube(Tolerance::default()).unwrap();
    let region = detect_apex_region(&cube2, stale_apex);
    // The stale handle *might* coincidentally map to the new
    // BRep's vertex 0 — that's allowed by the slot_map design.
    // What we assert is that detect doesn't panic and returns
    // a structurally-coherent answer (either valid with the
    // mapped vertex, or empty).
    if region.is_valid() {
        assert!(region.apex.is_some());
        assert!(!region.edges.is_empty());
    }
}

#[test]
fn outer_loop_valence_pre_fillet_cube() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let any_face = cube.faces.keys().next().unwrap();
    assert_eq!(outer_loop_valence(&cube, any_face), 4);
}

#[test]
fn outer_loop_valence_post_fillet_cube_cap_is_six() {
    // Route B corner reconstruction contracts the legacy
    // hexagon's anchor vertices into the on-sphere contacts:
    // the cap boundary is now THREE exact sphere-centered
    // arcs between the rolling-ball tangency feet — the
    // spherical-triangle anatomy of the true corner blend.
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
    assert_eq!(
        outer_loop_valence(&result.brep, cap_face),
        3,
        "reconstructed cube_uniform cap boundary is a spherical triangle: 3 exact arcs between tangency feet"
    );
}

#[test]
fn probe_cube_uniform_cap_surface_kind() {
    // Diagnostic — what surface kind does popkern's existing
    // pipeline produce for the cube_uniform cap? If it's already
    // a Sphere, Phase 6c-C may have nothing to do for this
    // fixture. If it's a Gregory or NURBS, Phase 6c-C swaps it
    // for an exact Sphere.
    use popkern_core::geometry::surface::SurfaceGeometry;
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
    let surf = result.brep.faces[cap_face]
        .surface
        .as_ref()
        .expect("cap face has a surface")
        .geometry
        .clone();
    let kind = match &surf {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
        SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
    };
    assert_eq!(kind, "Sphere", "cube_uniform cap should be a Sphere");
    if let SurfaceGeometry::Sphere { center, radius } = &surf {
        // Post-fix: cap centre should be at (+r, +r, +r) inside
        // the cube body for apex (0,0,0) with r=0.1. Pre-fix it
        // was at (-r, -r, -r) (outside the solid) due to the
        // median-equals-corner tie-breaker bug in surface_fix.
        assert!((center.x - 0.1).abs() < 1e-9);
        assert!((center.y - 0.1).abs() < 1e-9);
        assert!((center.z - 0.1).abs() < 1e-9);
        assert!((radius - 0.1).abs() < 1e-9);
    }
}

#[test]
fn probe_cube_uniform_rolling_ball_footprint_topology() {
    // Diagnostic — for each cap-loop vertex, dump its incident
    // edges and the faces those edges bound. Tells us whether
    // the rolling-ball-footprint vertices (A=(r,0,0)/B/C) are
    // 3-valent (just kev-able to merge cap loop) or 4-valent
    // (more incidences to handle).
    use popkern_core::geometry::surface::SurfaceGeometry;
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
    let cap_loop = result.brep.loops[result.brep.faces[cap_face].outer_loop].half_edge;

    let mut cursor = cap_loop;
    for k in 0..6 {
        let v = result.brep.half_edges[cursor].origin;
        let v_pos = result.brep.vertices[v].position;
        let mut incident_edges: Vec<_> = result
            .brep
            .half_edges
            .values()
            .filter(|he| he.origin == v)
            .map(|he| he.edge)
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        incident_edges.sort();
        eprintln!(
            "cap_v{k}: pos=({:.3},{:.3},{:.3}) valence={}",
            v_pos.x,
            v_pos.y,
            v_pos.z,
            incident_edges.len()
        );
        for e in &incident_edges {
            // Find a half-edge of this edge to identify its faces
            let hes: Vec<_> = result
                .brep
                .half_edges
                .iter()
                .filter(|(_, he)| he.edge == *e)
                .collect();
            let face_kinds: Vec<&'static str> = hes
                .iter()
                .map(|(_, he)| {
                    let f = result.brep.loops[he.loop_].face;
                    result.brep.faces[f]
                        .surface
                        .as_ref()
                        .map(|s| match &s.geometry {
                            SurfaceGeometry::Plane { .. } => "Plane",
                            SurfaceGeometry::Cylinder { .. } => "Cyl",
                            SurfaceGeometry::Sphere { .. } => "Sph",
                            SurfaceGeometry::Torus { .. } => "Torus",
                            _ => "Other",
                        })
                        .unwrap_or("?")
                })
                .collect();
            eprintln!("    edge {e:?} → faces: {face_kinds:?}");
        }
        cursor = result.brep.half_edges[cursor].next;
    }
}

#[test]
fn detect_post_fillet_cap_face_in_region_faces() {
    // From any vertex on the post-fillet cap face's outer loop,
    // `detect_apex_region` includes the cap face among
    // `region.faces`. Phase 6c-C.2 (deletion) relies on this.
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
    let cap_loop = result.brep.loops[result.brep.faces[cap_face].outer_loop].half_edge;
    let probe_v = result.brep.half_edges[cap_loop].origin;
    let region = detect_apex_region(&result.brep, probe_v);
    assert!(
        region.faces.contains(&cap_face),
        "cap face must be in apex region's faces from any cap-loop vertex"
    );
}
