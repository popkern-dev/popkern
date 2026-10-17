//! Trim-scale sphere and torus boxes contain their faces (PK-875): every
//! chart point the face's own region test calls inside lies in the box,
//! over primitives, boolean-trimmed caps holding either pole, and
//! rotated poses (pole caps off the world axes).

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Every inside sample of every sphere/torus face of `brep` lies in its
/// box; returns how many faces the trim-scale path boxed below their
/// whole surface.
fn check(brep: &BRep, what: &str) -> usize {
    let mut tighter = 0;
    for (fh, face) in &brep.faces {
        let Some(surface) = face.surface.as_ref().map(|s| &s.geometry) else {
            continue;
        };
        if !matches!(
            surface,
            SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Torus { .. }
        ) {
            continue;
        }
        let bx = crate::face_box::face_aabb(brep, fh).expect("a box");
        let ((u0, u1), (v0, v1)) = surface.parameter_range();
        let pad = 1e-7;
        let mut inside = 0;
        // Cell centres: a pole row is one 3-D point whose uv verdict
        // depends on an arbitrary u.
        for i in 0..97 {
            for j in 0..97 {
                let u = u0 + (u1 - u0) * (f64::from(i) + 0.5) / 97.0;
                let v = v0 + (v1 - v0) * (f64::from(j) + 0.5) / 97.0;
                let q = popkern_core::geometry::Point3::new(u, v, 0.0);
                if !crate::split_inner_loop::point_inside_face_region_uv(
                    brep,
                    fh,
                    q,
                    &brep.tolerance,
                ) {
                    continue;
                }
                inside += 1;
                let p = surface.point_at(u, v);
                assert!(
                    p.x >= bx.min.x - pad
                        && p.x <= bx.max.x + pad
                        && p.y >= bx.min.y - pad
                        && p.y <= bx.max.y + pad
                        && p.z >= bx.min.z - pad
                        && p.z <= bx.max.z + pad,
                    "{what}: face {fh:?} point {p:?} (uv {u:.3},{v:.3}) outside its box {bx:?}"
                );
            }
        }
        assert!(inside > 0, "{what}: face {fh:?} has no inside samples");
        let whole = match surface {
            SurfaceGeometry::Sphere { radius, .. } => 2.0 * radius,
            SurfaceGeometry::Torus {
                major_radius,
                minor_radius,
                ..
            } => 2.0 * (major_radius + minor_radius),
            _ => f64::INFINITY,
        };
        let size = (bx.max.x - bx.min.x)
            .max(bx.max.y - bx.min.y)
            .max(bx.max.z - bx.min.z);
        if size < 0.99 * whole {
            tighter += 1;
        }
    }
    tighter
}

fn posed(brep: BRep) -> BRep {
    let t = Transform3::translation(0.3, -0.2, 0.1)
        .compose(&Transform3::rotation_z(1.1))
        .compose(&Transform3::rotation_x(0.7));
    brep.transform(&t).expect("rigid")
}

#[test]
fn sphere_and_torus_face_boxes_contain_their_faces() {
    let sphere = BRep::sphere(1.0, tol()).expect("sphere");
    let torus = BRep::torus(1.0, 0.3, tol()).expect("torus");
    check(&sphere, "sphere");
    check(&torus, "torus");
    check(&posed(sphere.clone()), "posed sphere");
    check(&posed(torus.clone()), "posed torus");
    // Caps: the sphere cut by a slab keeps a cap at each pole.
    let slab = BRep::unit_cube(tol())
        .expect("cube")
        .transform(&Transform3::translation(-0.5, -0.5, -0.2))
        .expect("slab");
    let mut tighter = 0;
    for (what, body) in [
        (
            "sphere minus slab",
            crate::subtraction_ref(&sphere, &slab, tol(), None),
        ),
        (
            "sphere and slab",
            crate::intersection_ref(&sphere, &slab, tol(), None),
        ),
        (
            "torus minus slab",
            crate::subtraction_ref(&torus, &slab, tol(), None),
        ),
    ] {
        let body = body.unwrap_or_else(|e| panic!("{what}: {e}"));
        tighter += check(&body, what);
        tighter += check(&posed(body), what);
    }
    assert!(tighter > 0, "no face boxed below its whole surface");
}
