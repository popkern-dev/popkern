//! PK-327 tier-2 detection gates: the single-body self-intersection
//! scan finds transversal curves between non-adjacent faces and stays
//! silent on healthy bodies.

use popkern_bool::ops::self_intersect::find_self_intersections;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsSurface, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn healthy_cube_scans_clean() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let scan = find_self_intersections(&cube, &tol()).unwrap();
    assert!(scan.is_clean(), "clean cube must scan clean: {scan:?}");
    assert_eq!(scan.coincident_pairs_skipped, 0);
}

/// Unit cube whose top face is grafted with a biquadratic bulge dipping
/// through the bottom plane (center control at z = −5 ⇒ surface center
/// at z = −0.5). The boundary rows stay straight so the top face's
/// edges remain exactly on the grafted surface.
fn cube_with_piercing_top() -> BRep {
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let top = brep
        .faces
        .keys()
        .find(|&fh| {
            let ps = brep.face_outer_loop_positions(fh);
            !ps.is_empty() && ps.iter().all(|p| (p.z - 1.0).abs() < 1e-9)
        })
        .expect("top face");
    let z = |xi: usize, yi: usize| if xi == 1 && yi == 1 { -5.0 } else { 1.0 };
    let mut control_points = Vec::with_capacity(9);
    for yi in 0..3 {
        for xi in 0..3 {
            control_points.push(Point3::new(xi as f64 * 0.5, yi as f64 * 0.5, z(xi, yi)));
        }
    }
    brep.faces[top].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(NurbsSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points,
            num_u: 3,
            num_v: 3,
            weights: vec![1.0; 9],
        }),
        orientation: true,
    });
    brep.recompute_face_pcurves(top).expect("pcurves");
    brep
}

#[test]
fn piercing_bulge_is_detected_as_a_closed_interior_curve() {
    let brep = cube_with_piercing_top();
    let scan = find_self_intersections(&brep, &tol()).unwrap();
    assert!(
        !scan.curves.is_empty(),
        "bulge through the bottom plane must be detected (pairs checked: {})",
        scan.pairs_checked
    );
    // Every reported curve lies in the bottom plane (z = 0): the only
    // transversal pair is top-bulge × bottom.
    for c in &scan.curves {
        let (t0, t1) = c.curve.curve_3d.parameter_range();
        for i in 0..=8 {
            let t = t0 + (t1 - t0) * (i as f64) / 8.0;
            let p = c.curve.curve_3d.point_at(t);
            assert!(
                p.z.abs() < 1e-6,
                "scan curve leaves the bottom plane: {p:?}"
            );
            assert!(
                (0.0..=1.0).contains(&p.x) && (0.0..=1.0).contains(&p.y),
                "scan curve leaves the bottom face trim: {p:?}"
            );
        }
    }
}
