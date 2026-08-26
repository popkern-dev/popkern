//! End-to-end integration tests for popkern-shell.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Assert the Euler characteristic V - E + F = 2 for a simple solid (one shell,
/// genus 0, no inner loops beyond faces).
fn assert_euler_solid(brep: &BRep) {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    assert_eq!(
        v - e + f,
        2,
        "Euler failed for solid: V={v} E={e} F={f}, V-E+F={}",
        v - e + f
    );
}

/// Assert the Euler characteristic for a compound body (e.g., shell or thicken
/// result with outer + inner skins). The expected Euler characteristic is
/// V - E + F = 2 * C where C is the number of closed connected components.
/// For a shell/thicken with no removed faces, C = 2 (outer + inner).
fn assert_euler_compound(brep: &BRep, expected_components: i64) {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    let expected = 2 * expected_components;
    assert_eq!(
        v - e + f,
        expected,
        "Euler failed for compound body: V={v} E={e} F={f}, \
         V-E+F={} expected 2*{expected_components}={expected}",
        v - e + f
    );
}

/// Assert the body passes the FULL validator with zero errors. Served
/// shell results carry a structural serve gate in the pipeline (PK-321),
/// so tests hold them to the stronger validate_all bar; the historical
/// catch_unwind shim died with the stale-handle defect it tolerated.
fn assert_structurally_valid(brep: &BRep) {
    let errors = brep.validate_all();
    assert!(errors.is_empty(), "validation errors: {errors:?}");
}
#[test]
fn shell_cube_no_removal() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // Outer (6) + inner (6) = 12 faces, two disjoint closed surfaces.
    assert_eq!(result.brep.faces.len(), 12);
    assert_euler_compound(&result.brep, 2);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_cube_one_face_removed() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[first_face],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // Outer (5) + inner (5) + side walls (4) = 14 faces.
    assert_eq!(result.brep.faces.len(), 14);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_cube_two_adjacent_removed() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let faces: Vec<_> = brep.faces.keys().collect();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[faces[0], faces[1]],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // Outer (4) + inner (4) + side walls around the merged opening.
    // Two adjacent removed faces open one 6-edge rim: 4 + 4 + 6 = 14.
    assert_eq!(result.brep.faces.len(), 14);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_variable_thickness() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = shell_variable(
        &brep,
        0.1,
        &[(first_face, 0.2)],
        ShellDirection::Inward,
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    assert_eq!(result.brep.faces.len(), 12);
    assert_euler_compound(&result.brep, 2);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_inward_outward_midplane() {
    let brep = BRep::unit_cube(tol()).unwrap();

    let r_in = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    assert_eq!(r_in.brep.faces.len(), 12);
    assert_euler_compound(&r_in.brep, 2);
    assert_structurally_valid(&r_in.brep);

    let r_out = shell(
        &brep,
        0.1,
        ShellDirection::Outward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // Outward shell: the moved (void) skin is an outward offset —
    // rounded under the arc default (26 faces) + original skin (6).
    assert_eq!(r_out.brep.faces.len(), 32);
    assert_euler_compound(&r_out.brep, 2);
    assert_structurally_valid(&r_out.brep);

    // Midplane (PK-330): BOTH skins move — outer to ±0.05 outside the
    // unit cube, inner to 0.05..0.95 inside it. The pre-PK-321 pipeline
    // served an inward half-shell under this label.
    let r_mid = shell(
        &brep,
        0.1,
        ShellDirection::Midplane,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // Midplane: outward-moved skin rounds (26), inward skin stays
    // sharp (6) under the arc default.
    assert_eq!(r_mid.brep.faces.len(), 32);
    assert_euler_compound(&r_mid.brep, 2);
    assert_structurally_valid(&r_mid.brep);
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for (_, v) in &r_mid.brep.vertices {
        lo = lo.min(v.position.x);
        hi = hi.max(v.position.x);
    }
    assert!(
        (lo - (-0.05)).abs() < 1e-9 && (hi - 1.05).abs() < 1e-9,
        "outer midplane skin must sit at ±0.05 outside the cube: {lo}..{hi}"
    );
    let inner_xs: Vec<f64> = r_mid
        .brep
        .vertices
        .values()
        .map(|v| v.position.x)
        .filter(|&x| x > 0.0 && x < 1.0)
        .collect();
    assert_eq!(inner_xs.len(), 8, "8 inner-skin vertices strictly inside");
    for x in inner_xs {
        assert!(
            (x - 0.05).abs() < 1e-9 || (x - 0.95).abs() < 1e-9,
            "inner midplane skin at 0.05/0.95, got {x}"
        );
    }
}

#[test]
fn shell_midplane_with_removed_face_serves_rim_walls() {
    // Midplane with an opening: walls span the two MOVED skins and the
    // rim stays flush with the removed face's original plane.
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Midplane,
        &[face],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // PK-340: the outward-moved outer skin serves ROUNDED (8 gap
    // strips + 4 sphere corners at its closed corners + rim arcs at
    // the opening): outer 17 + inner 5 + wall quads 4 + rim corner
    // sectors 4 = 30 faces, one connected shell, no downgrades.
    assert_eq!(result.brep.faces.len(), 30);
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| matches!(w, popkern_shell::ShellWarning::JoinStyleDowngraded { .. })),
        "opening-rim corners must serve without downgrade: {:?}",
        result.warnings
    );
    assert_structurally_valid(&result.brep);
}

#[test]
fn offset_solid_outward() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result =
        popkern_shell::offset::offset_solid(&brep, 0.1, &ShellOptions::default(), tol()).unwrap();
    // Minkowski-faithful default (PK-325): outward offset of a convex
    // cube is the ROUNDED box — 6 planes + 12 cylinder strips + 8
    // sphere octant caps.
    assert_eq!(result.brep.faces.len(), 26);
    assert_eq!(result.brep.edges.len(), 48);
    assert_eq!(result.brep.vertices.len(), 24);
    assert_euler_solid(&result.brep);
    assert_structurally_valid(&result.brep);
    // Corner-sphere vertices sit at original corners displaced 0.1
    // along exactly one face normal: coordinates from {-0.1, 0, 1, 1.1}.
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                [-0.1, 0.0, 1.0, 1.1].iter().any(|e| (c - e).abs() < 1e-9),
                "vertex coord {c} unexpected after arc offset"
            );
        }
    }
}

#[test]
fn offset_solid_negative() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result =
        popkern_shell::offset::offset_solid(&brep, -0.1, &ShellOptions::default(), tol()).unwrap();
    // Topology counts must be preserved: 6 faces, 12 edges, 8 vertices.
    assert_eq!(result.brep.faces.len(), 6);
    assert_eq!(result.brep.edges.len(), 12);
    assert_eq!(result.brep.vertices.len(), 8);
    assert_euler_solid(&result.brep);
    assert_structurally_valid(&result.brep);
    // Every vertex must have moved from the original unit cube positions.
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                (c - 0.0).abs() > 1e-6 && (c - 1.0).abs() > 1e-6,
                "vertex coord {c} unchanged after negative offset"
            );
        }
    }
}

#[test]
fn thicken_cube() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result =
        popkern_shell::thicken::thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    // Rounded outer (26) + sharp inner (6) under the arc default
    // (no side walls for a closed solid).
    assert_eq!(result.brep.faces.len(), 32);
    assert_euler_compound(&result.brep, 2);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_empty_brep_error() {
    let brep = BRep::new(tol());
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    assert!(matches!(result, Err(ShellError::NoFaces)));
}

// ---- New tests below ----

#[test]
fn offset_solid_variable_per_face() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = popkern_shell::offset::offset_solid_variable(
        &brep,
        0.1,
        &[(first_face, 0.3)],
        &ShellOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(result.brep.faces.len(), 6);
    assert_euler_solid(&result.brep);
    assert_structurally_valid(&result.brep);
}

#[test]
fn thicken_prism_slab() {
    // Build a thin rectangular prism (acts as a slab) and thicken it.
    // The prism has 6 faces. Thicken produces outer + inner = 12 faces.
    use popkern_core::geometry::{Point3, Vec3};

    let profile = &[
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let brep = BRep::prism(profile, Vec3::new(0.0, 0.0, 1.0), 0.05, tol()).unwrap();
    // A 0.02 wall fits inside the 0.05-thick slab (0.01 inward per side).
    let result =
        popkern_shell::thicken::thicken(&brep, 0.02, &ShellOptions::default(), tol()).unwrap();
    // Thicken a closed solid: rounded outer (26) + sharp inner (6).
    assert_eq!(result.brep.faces.len(), 32);
    assert_euler_compound(&result.brep, 2);
    assert_structurally_valid(&result.brep);

    // A 0.2 wall does NOT fit: the inward half-offset (0.1) inverts the
    // 0.05-thick slab. This used to serve an inverted "hollow" body.
    let too_thick = popkern_shell::thicken::thicken(&brep, 0.2, &ShellOptions::default(), tol());
    assert!(
        matches!(too_thick, Err(ShellError::OffsetTooLarge { .. })),
        "expected OffsetTooLarge decline, got {too_thick:?}"
    );
}

#[test]
fn shell_cone_serves_with_exact_offset_arm() {
    // PK-322: the exact cone offset arm (same half-angle, apex
    // translated d/sin(half_angle) along the axis) plus the analytic
    // recompute arms serve a shelled cone validate_all-clean with
    // measured geometry — no NURBS approximation anywhere.
    use popkern_core::geometry::surface::SurfaceGeometry;

    let (radius, height, d) = (1.0_f64, 2.0_f64, 0.1_f64);
    let brep = BRep::cone(radius, height, tol()).unwrap();
    let base_cap = brep
        .faces
        .iter()
        .find(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            )
        })
        .map(|(h, _)| h)
        .unwrap();

    let result = shell(
        &brep,
        d,
        ShellDirection::Inward,
        &[base_cap],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();

    assert!(
        result.brep.validate_all().is_empty(),
        "shelled cone must be validate_all-clean: {:?}",
        result.brep.validate_all()
    );

    // The inner lateral face is a genuine Cone with the apex moved down
    // the axis by exactly d/sin(half_angle).
    let half_angle = (radius / height).atan();
    let expected_apex_z = height - d / half_angle.sin();
    let inner_apexes: Vec<f64> = result
        .brep
        .faces
        .values()
        .filter_map(|f| f.surface.as_ref())
        .filter_map(|s| match &s.geometry {
            SurfaceGeometry::Cone { apex, .. } => Some(apex.z),
            _ => None,
        })
        .filter(|z| (z - expected_apex_z).abs() < 1e-9)
        .collect();
    assert_eq!(
        inner_apexes.len(),
        1,
        "exactly one offset cone face with apex at z={expected_apex_z}"
    );

    // The inner apex VERTEX tracks the offset cone's apex exactly, and
    // the opening rim stays flush in the original base plane at the
    // inner cone's station radius.
    let apex_vertices = result
        .brep
        .vertices
        .values()
        .filter(|v| {
            v.position.x.abs() < 1e-9
                && v.position.y.abs() < 1e-9
                && (v.position.z - expected_apex_z).abs() < 1e-9
        })
        .count();
    assert_eq!(apex_vertices, 1, "inner apex vertex at the translated apex");

    let rim_radius = expected_apex_z * half_angle.tan();
    let rim_vertices = result
        .brep
        .vertices
        .values()
        .filter(|v| {
            let r = (v.position.x * v.position.x + v.position.y * v.position.y).sqrt();
            v.position.z.abs() < 1e-9 && (r - rim_radius).abs() < 1e-9
        })
        .count();
    assert!(
        rim_vertices >= 1,
        "opening rim vertex at radius {rim_radius} in the base plane"
    );
}
