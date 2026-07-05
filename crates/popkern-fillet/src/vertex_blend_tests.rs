//! Tests for §5.5 vertex blending via the public
//! [`vertex_blend_3d`] entry point. Exercises the
//! trihedral cube-corner case plus error paths.

use super::*;
use crate::types::FilletOptions;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Find the cube corner at `(0, 0, 0)` in a unit cube.
fn origin_corner(cube: &BRep) -> VertexHandle {
    cube.vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0 && v.position.z == 0.0)
        .map(|(h, _)| h)
        .unwrap()
}

fn count_face_types(brep: &BRep) -> (usize, usize, usize, usize) {
    let mut planes = 0;
    let mut fillets = 0;
    let mut spheres = 0;
    let mut others = 0;
    for fh in brep.faces.keys() {
        match brep.faces[fh].surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => planes += 1,
            Some(SurfaceGeometry::Torus { .. } | SurfaceGeometry::Cylinder { .. }) => fillets += 1,
            Some(SurfaceGeometry::Sphere { .. }) => spheres += 1,
            _ => others += 1,
        }
    }
    (planes, fillets, spheres, others)
}

/// Vertex blend at a single cube corner with radius 0.1.
/// Expect at least one sphere face appearing at the corner,
/// plus fillet strips on the three incident edges.
#[test]
fn vertex_blend_3d_on_cube_corner_produces_sphere_patch() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = origin_corner(&cube);
    let result = vertex_blend_3d(&cube, corner, 0.1, &FilletOptions::default(), tol(), None)
        .expect("vertex blend on cube corner should succeed");

    let (planes, fillets, spheres, _others) = count_face_types(&result.brep);
    assert_eq!(
        fillets, 3,
        "expected 3 fillet strips (one per incident edge), got {fillets}"
    );
    assert_eq!(
        spheres, 1,
        "expected 1 sphere patch at the blended corner, got {spheres}"
    );
    assert!(
        planes >= 6,
        "expected ≥6 planar faces (3 trimmed cube faces + 3 untouched on the far side), got {planes}"
    );
    result.brep.validate().unwrap();
}

/// The sphere patch's radius must equal the requested blend
/// radius (within tolerance) — the pipeline's corner-resolution
/// uses a rolling-ball tangent sphere at exactly that radius.
#[test]
fn vertex_blend_3d_sphere_radius_matches_request() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = origin_corner(&cube);
    let r = 0.1_f64;
    let result = vertex_blend_3d(&cube, corner, r, &FilletOptions::default(), tol(), None).unwrap();

    let sphere_radius = result
        .brep
        .faces
        .values()
        .find_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Sphere { radius, .. }) => Some(*radius),
            _ => None,
        })
        .expect("no sphere face found");
    assert!(
        (sphere_radius - r).abs() < 1e-6,
        "sphere radius {sphere_radius} != requested {r}"
    );
}

/// Invalid vertex handle is rejected with `InvalidParameter`.
#[test]
fn vertex_blend_3d_rejects_invalid_vertex_handle() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let bogus = VertexHandle::default();
    let err =
        vertex_blend_3d(&cube, bogus, 0.1, &FilletOptions::default(), tol(), None).unwrap_err();
    assert!(matches!(
        err,
        BlendError::InvalidParameter { ref parameter, .. } if parameter == "vertex"
    ));
}

/// Non-finite radius is rejected at the API boundary.
#[test]
fn vertex_blend_3d_rejects_invalid_radius() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = origin_corner(&cube);
    for bad in [-0.1, 0.0, f64::NAN, f64::INFINITY] {
        let err = vertex_blend_3d(&cube, corner, bad, &FilletOptions::default(), tol(), None)
            .unwrap_err();
        assert!(
            matches!(err, BlendError::InvalidParameter { ref parameter, .. } if parameter == "radius"),
            "radius={bad} should produce InvalidParameter(radius)"
        );
    }
}

/// Vertex blend on the same corner at two different radii
/// produces sphere faces with correspondingly different radii.
#[test]
fn vertex_blend_3d_radius_scales_sphere_patch() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = origin_corner(&cube);
    let small =
        vertex_blend_3d(&cube, corner, 0.05, &FilletOptions::default(), tol(), None).unwrap();
    let big = vertex_blend_3d(&cube, corner, 0.15, &FilletOptions::default(), tol(), None).unwrap();

    let small_r = small
        .brep
        .faces
        .values()
        .find_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Sphere { radius, .. }) => Some(*radius),
            _ => None,
        })
        .unwrap();
    let big_r = big
        .brep
        .faces
        .values()
        .find_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Sphere { radius, .. }) => Some(*radius),
            _ => None,
        })
        .unwrap();
    assert!((small_r - 0.05).abs() < 1e-6);
    assert!((big_r - 0.15).abs() < 1e-6);
}

// ---------------------------------------------------------------------------
// Follow-up: G1 continuity at the blend boundary
// ---------------------------------------------------------------------------

/// Rolling-ball invariant for a cube corner: the sphere's
/// centre must sit at distance `r · √3` from the corner along
/// the interior diagonal, and its outward normals at the three
/// tangent points `(0, r, r)`, `(r, 0, r)`, `(r, r, 0)` must
/// match the three cube faces' outward normals (planes `x = 0`,
/// `y = 0`, `z = 0`, all pointing in −x / −y / −z). Checking
/// the centre position gives the tightest geometric invariant
/// available from the final surface data — a rolling-ball
/// blend with radius r at a trihedral right-angle corner has
/// exactly one valid sphere centre, and that centre is
/// `(r, r, r)` for the origin corner of a unit cube.
#[test]
fn vertex_blend_3d_sphere_centre_at_rolling_ball_position() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = origin_corner(&cube);
    let r = 0.1_f64;
    let result = vertex_blend_3d(&cube, corner, r, &FilletOptions::default(), tol(), None).unwrap();

    let brep = &result.brep;
    let (sphere_centre, sphere_radius) = brep
        .faces
        .values()
        .find_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Sphere { center, radius }) => Some((*center, *radius)),
            _ => None,
        })
        .expect("no sphere face found");

    assert!((sphere_radius - r).abs() < 1e-9);

    // Distance from the corner to the sphere centre must be
    // exactly `r · √3` (the rolling-ball invariant for a
    // right-angle trihedral corner), and the offset must be
    // along the ±(1, 1, 1) body diagonal (each component has
    // absolute value `r`). The sign depends on which side of
    // the corner the cube's interior lies on — `unit_cube`'s
    // interior may extend in either direction from (0, 0, 0),
    // so we don't prescribe it.
    let corner_pos = cube.vertices[corner].position;
    let offset = [
        sphere_centre.x - corner_pos.x,
        sphere_centre.y - corner_pos.y,
        sphere_centre.z - corner_pos.z,
    ];
    let diag = (offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2]).sqrt();
    let expected_diag = r * 3.0_f64.sqrt();
    assert!(
        (diag - expected_diag).abs() < 1e-9,
        "corner-to-sphere-centre distance {diag} != expected r·√3 = {expected_diag}"
    );
    for c in offset {
        assert!(
            (c.abs() - r).abs() < 1e-9,
            "sphere-centre offset component {c} should have |value| = r = {r}"
        );
    }
}

// ---------------------------------------------------------------------------
// Follow-up: Non-trihedral vertex (4-valence pyramid apex)
// ---------------------------------------------------------------------------

/// A square pyramid's apex has 4 incident edges (one per
/// slant), so vertex-blending it exercises the N-gon path in
/// `insert_vertex_blend` that didn't get covered by the cube
/// tests (all 3-valence). Expect 4 fillet strips plus per-fillet
/// sphere patches at the apex (Phase 6c-E: per-fillet rolling-ball
/// dispatcher places one sphere centred on each cylinder axis).
#[test]
#[ignore = "mirror-era baseline: pre-PK-099 ball placement put pyramid strips on the exterior mirror; the apex Direct/Setback/Gregory rungs were derived against that geometry and need re-derivation on real inscribed strips (PK-101, root cause of the PK-092 direct_apex 8/8 baseline)"]
fn vertex_blend_3d_on_pyramid_apex_produces_sphere_patch() {
    let pyramid = BRep::square_pyramid(2.0, 1.5, tol()).unwrap();
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0 && v.position.z > 0.0)
        .map(|(h, _)| h)
        .unwrap();

    let result = vertex_blend_3d(&pyramid, apex, 0.1, &FilletOptions::default(), tol(), None)
        .expect("vertex blend on pyramid apex (4-valence) should succeed");

    let (_planes, tori, spheres, _others) = count_face_types(&result.brep);
    assert_eq!(
        tori, 4,
        "expected 4 fillet strips (one per slant edge), got {tori}"
    );
    assert!(
        spheres >= 1,
        "expected at least one sphere patch at the apex, got {spheres}"
    );
    result.brep.validate().unwrap();
}

// ---------------------------------------------------------------------------
// Follow-up: Curved-face support
// ---------------------------------------------------------------------------

/// Vertex blending works when adjacent faces are curved. The
/// cleanest properly-constructed test vehicle is the 12-edge
/// cube fillet: every sphere vertex-blend face is bounded by
/// **fillet strips** (genuinely curved surfaces), not planes. If
/// the batch pipeline's corner-resolution couldn't handle
/// curved adjacencies, the torus-sphere boundaries would fail.
///
/// This test deliberately goes through the same `fillet_edges`
/// path that `vertex_blend_3d` delegates to, and verifies that
/// all 8 cube corners produce sphere patches bounded by tori.
#[test]
fn vertex_blend_3d_curved_face_support_via_full_cube_fillet() {
    use crate::fillet_edges;
    use crate::types::{FilletRadius, FilletSpec};
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    let (_, tori, spheres, _) = count_face_types(&result.brep);
    assert_eq!(tori, 12, "expected 12 torus strip faces");
    assert_eq!(
        spheres, 8,
        "each cube corner should have a sphere vertex-blend face bounded by fillet strips"
    );
    result.brep.validate().unwrap();
}

// ---------------------------------------------------------------------------
// Follow-up: Variable-radius per edge
// ---------------------------------------------------------------------------

/// Variable radii across the three edges at a cube corner:
/// 0.05, 0.10, 0.15. The batch pipeline should produce 3 torus
/// strips with correspondingly different radii and a sphere
/// patch at the corner whose radius is driven by the local
/// rolling-ball geometry (not necessarily equal to any single
/// edge radius, but a function of the mixed spec).
#[test]
fn vertex_blend_3d_variable_radii_produce_distinct_torus_sizes() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = origin_corner(&cube);
    // Collect the three edges incident to the origin corner.
    let incident: Vec<_> = cube
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == corner)
        .map(|(_, he)| he.edge)
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();
    assert_eq!(incident.len(), 3);

    let edge_radii: Vec<_> = incident
        .iter()
        .copied()
        .zip([0.05_f64, 0.10, 0.15])
        .collect();
    let result = vertex_blend_3d_variable(
        &cube,
        corner,
        &edge_radii,
        &FilletOptions::default(),
        tol(),
        None,
    )
    .expect("variable-radius vertex blend should succeed");

    // Verify we got three tori with distinct minor radii
    // matching the spec. The sphere vertex blend is OPTIONAL
    // here: a rolling-ball sphere tangent to all three tori
    // only exists when the three radii are equal, so mixed
    // radii may either (a) skip the sphere patch, (b) pick a
    // single compromise sphere radius, or (c) fail. We accept
    // all three outcomes as long as the torus geometry is
    // correct and the BRep still validates.
    let (_, tori, _spheres, _) = count_face_types(&result.brep);
    assert_eq!(tori, 3, "expected 3 fillet strips, got {tori}");

    let mut minor_radii: Vec<f64> = result
        .brep
        .faces
        .values()
        .filter_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Torus { minor_radius, .. }) => Some(*minor_radius),
            Some(SurfaceGeometry::Cylinder { radius, .. }) => Some(*radius),
            _ => None,
        })
        .collect();
    minor_radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(minor_radii.len(), 3);
    let expected = [0.05, 0.10, 0.15];
    for (got, want) in minor_radii.iter().zip(expected.iter()) {
        assert!(
            (got - want).abs() < 1e-6,
            "torus minor_radius {got} != {want}"
        );
    }
    result.brep.validate().unwrap();
}

/// Variable-radius API rejects non-finite radii.
#[test]
fn vertex_blend_3d_variable_rejects_bad_radii() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = origin_corner(&cube);
    let first_edge = cube.edges.keys().next().unwrap();
    let err = vertex_blend_3d_variable(
        &cube,
        corner,
        &[(first_edge, f64::NAN)],
        &FilletOptions::default(),
        tol(),
        None,
    )
    .unwrap_err();
    assert!(matches!(
        err,
        BlendError::InvalidParameter { ref parameter, .. } if parameter == "edge_radii"
    ));
}

/// Variable-radius API rejects an edge that isn't incident to
/// the target vertex.
#[test]
fn vertex_blend_3d_variable_rejects_non_incident_edge() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = origin_corner(&cube);
    // Find an edge NOT incident to `corner`. A unit cube has
    // 12 edges; 3 touch each corner, so 9 don't touch `corner`.
    let far_edge = cube
        .edges
        .iter()
        .find(|(_, e)| {
            let he = e.half_edge;
            let a = cube.half_edges[he].origin;
            let b = cube.half_edges[cube.half_edges[he].twin].origin;
            a != corner && b != corner
        })
        .map(|(h, _)| h)
        .unwrap();
    let err = vertex_blend_3d_variable(
        &cube,
        corner,
        &[(far_edge, 0.1)],
        &FilletOptions::default(),
        tol(),
        None,
    )
    .unwrap_err();
    assert!(matches!(
        err,
        BlendError::InvalidParameter { ref parameter, .. } if parameter == "edge_radii"
    ));
}

/// Variable-radius API rejects a partial spec that omits some
/// incident edges of `vertex`.
#[test]
fn vertex_blend_3d_variable_rejects_partial_incident_spec() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = origin_corner(&cube);
    let one_edge = cube
        .half_edges
        .iter()
        .find(|(_, he)| he.origin == corner)
        .map(|(_, he)| he.edge)
        .unwrap();
    // Only one of the three incident edges is specified.
    let err = vertex_blend_3d_variable(
        &cube,
        corner,
        &[(one_edge, 0.1)],
        &FilletOptions::default(),
        tol(),
        None,
    )
    .unwrap_err();
    assert!(matches!(
        err,
        BlendError::InvalidParameter { ref parameter, .. } if parameter == "edge_radii"
    ));
}
