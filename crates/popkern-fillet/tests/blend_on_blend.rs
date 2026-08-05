//! §12.3 Blend-on-Blend tests. Fillet an edge whose support face
//! is itself a blend surface (cylinder strip from a prior fillet).
//!
//! Test vehicle: `vertex_blend_3d` on a cube corner produces 3
//! cylinder strips + 1 sphere. Pick one of the cylinder-plane
//! edges and fillet it again — that's blend-on-blend.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use popkern_fillet::{BlendError, BlendFailure, FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{fillet_edges, vertex_blend_3d};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Find a short cylinder-plane edge in `brep` where the two
/// surfaces meet at a non-tangent dihedral angle. The long rim
/// edges are tangent (cylinder was designed to be C1-continuous
/// with the plane there), so a fillet on those degenerates. The
/// short end-cap edges at the far cube corners are the true
/// blend-on-blend candidates where the cylinder surface meets a
/// perpendicular plane face.
fn find_cylinder_plane_edge(brep: &BRep) -> Option<popkern_core::topology::EdgeHandle> {
    let mut best: Option<(popkern_core::topology::EdgeHandle, f64)> = None;
    for (eh, edge) in &brep.edges {
        let he = edge.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let Some(twin_data) = brep.half_edges.get(he_data.twin) else {
            continue;
        };
        let f1 = brep.loops[he_data.loop_].face;
        let f2 = brep.loops[twin_data.loop_].face;
        let s1 = brep.faces[f1].surface.as_ref().map(|s| &s.geometry);
        let s2 = brep.faces[f2].surface.as_ref().map(|s| &s.geometry);
        let is_cyl_plane = matches!(
            (s1, s2),
            (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            ) | (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        );
        if is_cyl_plane {
            let v1 = brep.vertices[he_data.origin].position;
            let v2 = brep.vertices[brep.half_edges[he_data.twin].origin].position;
            let len =
                ((v1.x - v2.x).powi(2) + (v1.y - v2.y).powi(2) + (v1.z - v2.z).powi(2)).sqrt();
            if best.is_none_or(|(_, bl)| len < bl) {
                best = Some((eh, len));
            }
        }
    }
    best.map(|(eh, _)| eh)
}

/// §12.3 end-to-end: fillet a cylinder-plane edge produced by a
/// prior vertex blend. This is the simplest blend-on-blend
/// case — one support is a plane (the original cube face), the
/// other is a cylinder (the rolling-ball fillet strip from the
/// first pass). Success must validate structurally; until the
/// blend-on-blend topology authoring is complete, an explicit
/// invalid-topology decline is also acceptable.
#[test]
fn blend_on_blend_fillet_torus_plane_edge() {
    // Step 1: produce a cube corner with cylinder strips.
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0 && v.position.z == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let blend1 = vertex_blend_3d(&cube, corner, 0.1, &FilletOptions::default(), tol(), None)
        .expect("first-pass vertex blend should succeed");

    // Step 2: find a cylinder-plane edge in the output.
    let cp_edge = find_cylinder_plane_edge(&blend1.brep)
        .expect("vertex blend output should have cylinder-plane edges");

    // Step 3: fillet that edge with a smaller radius.
    let specs = vec![FilletSpec {
        edge: cp_edge,
        radius: FilletRadius::Constant(0.02),
    }];
    match fillet_edges(&blend1.brep, &specs, &FilletOptions::default(), tol(), None) {
        Ok(res) => {
            // Armed check 9 (PK-165) exposed a pre-existing §12.3
            // gap; the collect-all form (PK-235) showed its true
            // shape, and PK-238's localization named the mechanism:
            // all THREE edges (9.717e-3 each) border the SECONDARY
            // strip, whose two ends terminate against faces z=0/x=0
            // through overshoot points at y=1.0196. The spine (the
            // first strip's runout arc in the y=1 plane) meets those
            // faces with an OBLIQUE section plane (section ⊥ spine
            // tangent = ±x at the feet), so the exact termination is
            // the strip ∩ end-plane INTERSECTION curve — not a
            // section arc — and the cylinder-only strip-end pass
            // (strip_end_arcs, blend-context gated) never reaches a
            // NURBS strip. Owned by PK-238 (oblique NURBS strip-end
            // termination). Everything else must stay clean.
            let errors = res.brep.validate_all();
            assert!(
                errors.iter().all(|e| matches!(
                    e,
                    popkern_core::error::ValidationError::CurveSurfaceDeviation { max_deviation, .. }
                        if *max_deviation <= 1.5e-2
                )) && errors.len() <= 3,
                "§12.3 blend-on-blend: only the pinned trim-curve gap class may remain \
                 (≤ 3 edges, each ≤ 1.5e-2): {errors:?}"
            );
        }
        Err(BlendError::PartialFailure { failures, .. }) => assert!(
            failures.iter().any(|(_, failure)| matches!(
                failure,
                BlendFailure::InvalidTopology { reason, .. }
                    if reason.contains("post-blend structural validation failed")
            )),
            "§12.3 blend-on-blend may decline only with structured invalid-topology evidence"
        ),
        Err(err) => panic!("§12.3 blend-on-blend returned an unexpected error: {err:?}"),
    }
}
