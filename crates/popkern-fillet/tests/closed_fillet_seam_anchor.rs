//! Closed-fillet seam-anchor consistency pin (PK-285).
//!
//! A closed rim fillet on a seamed support (cylinder side face) must
//! move the rim-seam vertex to a SEAM-CONSISTENT point on the contact
//! circle — the point radially nearest the old vertex — never to the
//! trim curve's own parameterization origin. That origin is an
//! offset-surface-intersection marching artifact: ULP-level input
//! changes (measured: a 7e-15 extrude-height difference from a
//! derivative-evaluation change) moved it to the diametrically
//! opposite side of the cylinder, and the rebuilt seam edge became a
//! chord crossing the solid — caught by check 9 at deviation ≈ r.
//!
//! This pin builds the boss-on-plate class at kernel level (block ∪
//! cylinder), fillets BOTH Cylinder×Plane rims (base seam and top
//! rim), and asserts the result is check-9 clean and every
//! cylinder-seam edge endpoint still lies on its cylinder.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[test]
fn two_offset_bosses_fillet_watertight() {
    // The double-boss class: two closed base-seam fillets in ONE
    // operation, boss centers OFFSET from the support's chart origin
    // (the corpus double_boss_plate regression class).
    let tol = Tolerance::default();
    let plate = BRep::unit_cube(tol).expect("cube");
    let plate = plate
        .transform(&Transform3::scale_uniform(80.0))
        .expect("scale");
    let boss_a = BRep::cylinder(8.0, 12.0, tol)
        .expect("boss a")
        .transform(&Transform3::translation(20.0, 40.0, 80.0))
        .expect("translate a");
    let boss_b = BRep::cylinder(6.0, 15.0, tol)
        .expect("boss b")
        .transform(&Transform3::translation(60.0, 40.0, 80.0))
        .expect("translate b");
    let body = boolean(BoolOperation::Union, plate, boss_a, tol, None).expect("union a");
    let body = boolean(BoolOperation::Union, body, boss_b, tol, None).expect("union b");

    let specs: Vec<FilletSpec> = body
        .edges
        .keys()
        .filter(|&eh| {
            let e = &body.edges[eh];
            let he = &body.half_edges[e.half_edge];
            let f1 = body.loops[he.loop_].face;
            let f2 = body.loops[body.half_edges[he.twin].loop_].face;
            let s1 = body.faces[f1].surface.as_ref().map(|s| &s.geometry);
            let s2 = body.faces[f2].surface.as_ref().map(|s| &s.geometry);
            // Base seams only: cylinder × the plate top plane.
            let is_pair = matches!(
                (s1, s2),
                (
                    Some(SurfaceGeometry::Cylinder { .. }),
                    Some(SurfaceGeometry::Plane { .. })
                ) | (
                    Some(SurfaceGeometry::Plane { .. }),
                    Some(SurfaceGeometry::Cylinder { .. })
                )
            );
            let he2 = &body.half_edges[e.half_edge];
            let z_ok = body.vertices[he2.origin].position.z < 80.5;
            is_pair && z_ok
        })
        .map(|eh| FilletSpec {
            edge: eh,
            radius: FilletRadius::Constant(2.0),
        })
        .collect();
    assert_eq!(specs.len(), 2, "two base seams expected");

    let result = fillet_edges(&body, &specs, &FilletOptions::default(), tol, None)
        .expect("double boss fillets succeed");
    let deviations = result.brep.validate_curve_surface_all();
    assert!(
        deviations.is_empty(),
        "curve-surface deviations: {deviations:?}"
    );
    result
        .brep
        .validate()
        .expect("double boss fillet result validates");
}

#[test]
fn closed_rim_fillets_keep_seam_edges_on_surface() {
    let tol = Tolerance::default();
    let block = BRep::unit_cube(tol).expect("cube");
    let block = block
        .transform(&Transform3::scale_uniform(50.0))
        .expect("scale");
    let boss = BRep::cylinder(10.0, 15.0, tol).expect("cylinder");
    let boss = boss
        .transform(&Transform3::translation(25.0, 25.0, 50.0))
        .expect("translate");
    let body = boolean(BoolOperation::Union, block, boss, tol, None).expect("union");

    let specs: Vec<FilletSpec> = body
        .edges
        .keys()
        .filter(|&eh| {
            let e = &body.edges[eh];
            let he = &body.half_edges[e.half_edge];
            let f1 = body.loops[he.loop_].face;
            let f2 = body.loops[body.half_edges[he.twin].loop_].face;
            let s1 = body.faces[f1].surface.as_ref().map(|s| &s.geometry);
            let s2 = body.faces[f2].surface.as_ref().map(|s| &s.geometry);
            matches!(
                (s1, s2),
                (
                    Some(SurfaceGeometry::Cylinder { .. }),
                    Some(SurfaceGeometry::Plane { .. })
                ) | (
                    Some(SurfaceGeometry::Plane { .. }),
                    Some(SurfaceGeometry::Cylinder { .. })
                )
            )
        })
        .map(|eh| FilletSpec {
            edge: eh,
            radius: FilletRadius::Constant(3.0),
        })
        .collect();
    assert_eq!(specs.len(), 2, "base seam + top rim expected");

    let result = fillet_edges(&body, &specs, &FilletOptions::default(), tol, None)
        .expect("closed rim fillets succeed");

    // Check 9 across the whole result: no edge curve off its face
    // surfaces (the chord regression reported deviation ≈ r here).
    let deviations = result.brep.validate_curve_surface_all();
    assert!(
        deviations.is_empty(),
        "curve-surface deviations after closed rim fillets: {deviations:?}"
    );

    // Every seam edge of a cylinder face (face_1 == face_2) must have
    // both endpoints ON that cylinder — the direct signature of the
    // wrong-anchor class.
    for (eh, e) in &result.brep.edges {
        let he = &result.brep.half_edges[e.half_edge];
        let twin = &result.brep.half_edges[he.twin];
        let f1 = result.brep.loops[he.loop_].face;
        let f2 = result.brep.loops[twin.loop_].face;
        if f1 != f2 {
            continue;
        }
        let Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) = result.brep.faces[f1].surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        for v in [he.origin, twin.origin] {
            let p = result.brep.vertices[v].position;
            let d = p - *origin;
            let radial = d - *axis * d.dot(axis);
            assert!(
                (radial.length() - radius).abs() < 1e-6,
                "seam edge {eh:?} endpoint {p:?} off its cylinder (radial {} vs r {radius})",
                radial.length()
            );
        }
    }
}
