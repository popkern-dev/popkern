//! Vanishing fillet ends (PK-764, `docs/architecture/fillet-vanishing-ends.md`):
//! the equal-radius axis tee's saddle edges end at the singular points,
//! where the cylinders are tangent. Each saddle strip computes. Its walk
//! starts at the regular end, both rails end exactly on the tangent
//! vertex, and the strip chart stops there at zero width instead of
//! folding past it.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// The axis tee: an r 40 main cylinder about z, an r 40 branch along
/// +x from the main axis at z = 60.
fn axis_tee(tol: Tolerance) -> BRep {
    let main = BRep::cylinder(40.0, 120.0, tol).expect("main");
    let branch = BRep::cylinder(40.0, 100.0, tol)
        .expect("branch")
        .transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2))
        .expect("turn")
        .transform(&Transform3::translation(0.0, 0.0, 60.0))
        .expect("place");
    popkern_bool::union(main, branch, tol, None).expect("union")
}

#[test]
fn equal_radius_saddle_strips_run_out_at_the_tangent_points() {
    let tol = Tolerance::default();
    let tee = axis_tee(tol);
    let singular = [Point3::new(0.0, 40.0, 60.0), Point3::new(0.0, -40.0, 60.0)];
    let mut saddle = 0;
    for (edge, _) in tee.edges.iter() {
        let Ok(mut analyses) = crate::analysis::analyze_edges(&tee, &[edge], false) else {
            continue;
        };
        let a = analyses.remove(0);
        if !(matches!(a.surface_1, SurfaceGeometry::Cylinder { .. })
            && matches!(a.surface_2, SurfaceGeometry::Cylinder { .. })
            && a.face_1 != a.face_2)
        {
            continue;
        }
        saddle += 1;
        let g = crate::geometry_compute::compute_fillet_geometry(
            &a,
            &crate::types::FilletRadius::Constant(4.0),
            &tee,
            &tol,
            crate::types::ContinuityLevel::G1,
            crate::types::CapKind::Open,
            crate::types::BlendProfile::default(),
        )
        .unwrap_or_else(|e| panic!("saddle edge {edge:?}: {e:?}"));
        // Both rails end on the edge's tangent vertex.
        let ends = |c: &popkern_core::geometry::NurbsCurve| {
            let (s, e) = c.parameter_range();
            [c.point_at(s), c.point_at(e)]
        };
        for rail in [&g.trim_3d_1, &g.trim_3d_2] {
            let at_tangent = ends(rail)
                .iter()
                .any(|p| singular.iter().any(|s| p.distance_to(s) < 1e-6));
            assert!(at_tangent, "edge {edge:?}: a rail misses the tangent point");
        }
        // The chart stops at the tangent point: its end row is the
        // point itself, and nothing lies past it (x < 0 is beyond).
        let ((u0, u1), (v0, v1)) = g.surface.parameter_range();
        let end_row = [u0, u1].into_iter().find(|&u| {
            let p = g.surface.point_at(u, v0);
            singular.iter().any(|s| p.distance_to(s) < 1e-6)
        });
        let u = end_row.unwrap_or_else(|| panic!("edge {edge:?}: no collapsed end row"));
        let width = g
            .surface
            .point_at(u, v0)
            .distance_to(&g.surface.point_at(u, v1));
        assert!(width < 1e-6, "edge {edge:?}: end row width {width}");
        for k in 0..=20 {
            let p = g
                .surface
                .point_at(u0 + (u1 - u0) * k as f64 / 20.0, 0.5 * (v0 + v1));
            assert!(
                p.x > -1e-6,
                "edge {edge:?}: the strip runs past the tangent point: {p:?}"
            );
        }
    }
    assert_eq!(saddle, 4, "the axis tee has four saddle edges");
}
