//! Planar regression and pcurve preservation tests for boolean operations
//! involving curved surfaces.
//!
//! Split from integration_tests_curved.rs to stay under the 500-line limit.

// Structural lint marker: this file is test code, included via #[cfg(test)] in lib.rs.
#[cfg(test)]
use popkern_core::{
    geometry::{Transform3, Vec3},
    tolerance::Tolerance,
    topology::BRep,
};

/// Translate all vertices and surface/curve geometry of a BRep by `offset`.
fn translate(mut brep: BRep, offset: Vec3) -> BRep {
    use popkern_core::geometry::{curve::CurveGeometry, surface::SurfaceGeometry};
    for (_, vertex) in &mut brep.vertices {
        vertex.position = vertex.position + offset;
    }
    for (_, face) in &mut brep.faces {
        if let Some(ref mut sof) = face.surface {
            match sof.geometry {
                SurfaceGeometry::Plane { ref mut origin, .. } => *origin = *origin + offset,
                SurfaceGeometry::Cylinder { ref mut origin, .. } => *origin = *origin + offset,
                SurfaceGeometry::Cone { ref mut apex, .. } => *apex = *apex + offset,
                SurfaceGeometry::Sphere { ref mut center, .. } => *center = *center + offset,
                _ => {}
            }
        }
    }
    for (_, edge) in &mut brep.edges {
        if let Some(ref mut coe) = edge.curve {
            match coe.geometry {
                CurveGeometry::Line { ref mut origin, .. } => *origin = *origin + offset,
                CurveGeometry::Circle { ref mut center, .. }
                | CurveGeometry::Arc { ref mut center, .. } => *center = *center + offset,
                _ => {}
            }
        }
    }
    brep
}

// ── Planar regression: cube-cube booleans still validate after curved changes ─

#[test]
fn cube_cube_union_still_validates() {
    let t = Tolerance::default();
    let a = BRep::unit_cube(t.clone()).unwrap();
    let b = translate(
        BRep::unit_cube(t.clone()).unwrap(),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let result = crate::union(a, b, t, None).unwrap();
    let errs = result.validate_all();
    assert!(
        errs.is_empty(),
        "Cube-cube union should still validate: {:?}",
        errs
    );
}

#[test]
fn overlapping_cube_subtraction_still_validates() {
    let t = Tolerance::default();
    let a = BRep::unit_cube(t.clone()).unwrap();
    let b = translate(
        BRep::unit_cube(t.clone()).unwrap(),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let result = crate::subtraction(a, b, t, None).unwrap();
    let errs = result.validate_all();
    assert!(
        errs.is_empty(),
        "Cube-cube subtraction should validate: {:?}",
        errs
    );
}

// ── Pcurve preservation through boolean pipeline ─────────────────────────────

#[test]
fn disjoint_cube_union_preserves_pcurves() {
    let t = Tolerance::default();
    let a = BRep::unit_cube(t.clone()).unwrap();
    let b = translate(
        BRep::unit_cube(t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::union(a, b, t, None).unwrap();
    let with_pcurves = result
        .half_edges
        .values()
        .filter(|he| he.pcurve.is_some())
        .count();
    assert!(
        with_pcurves > 0,
        "Result should have pcurves from cube constructors"
    );
}

#[test]
fn boolean_subtraction_preserves_and_creates_pcurves() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.2, 2.0, tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let result = crate::subtraction(cube, cyl, tol, None).unwrap();

    let total = result.half_edges.len();
    let with_pcurves = result
        .half_edges
        .values()
        .filter(|he| he.pcurve.is_some())
        .count();

    // Most half-edges should have pcurves (cube edges from operand + intersection edges).
    assert!(
        with_pcurves > total / 2,
        "expected most half-edges to have pcurves, got {with_pcurves}/{total}"
    );
}
