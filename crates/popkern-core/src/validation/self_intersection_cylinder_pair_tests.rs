//! PK-866 regressions for the cylinder × cylinder fast path: a real
//! crossing of two trimmed cylinder walls is flagged, the same walls
//! moved apart are not, and near samples refine onto the crossing.

use crate::error::ValidationError;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};

/// A four-sided face through `corners` carrying `surface`.
fn quad_face(brep: &mut BRep, corners: [Point3; 4], surface: SurfaceGeometry) -> FaceHandle {
    let mvfs = brep.mvfs(corners[0]).unwrap();
    let e1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, corners[1])
        .unwrap();
    let e2 = brep.mev(e1.half_edge_twin, corners[2]).unwrap();
    let e3 = brep.mev(e2.half_edge_twin, corners[3]).unwrap();
    let face = brep.mef(e3.half_edge_twin, e1.half_edge).unwrap().face;
    brep.faces[face].surface = Some(SurfaceOnFace {
        geometry: surface,
        orientation: true,
    });
    face
}

/// Quarter wall of radius 0.5 about +Z (z in [0, 1]), and a wall of
/// radius 0.4 about +X through (0, 0, 0.5) spanning x in `x_range`,
/// angle in [-pi/4, pi/4] about its axis.
fn two_walls(x_range: (f64, f64)) -> (BRep, FaceHandle, FaceHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let a = quad_face(
        &mut brep,
        [
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
            Point3::new(0.0, 0.5, 1.0),
            Point3::new(0.5, 0.0, 1.0),
        ],
        SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 0.5,
        },
    );
    // Corners from B's own chart, just off the +Y direction (its u
    // seam sits there), so the positional trim window spans the
    // quarter rather than wrapping the long way round.
    use crate::geometry::surface::SurfaceEval;
    let surface_b = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.5),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 0.4,
    };
    let (u_mid, _) = surface_b
        .closest_point(&Point3::new(0.3, 0.4, 0.5))
        .unwrap();
    let (u0, u1) = (u_mid + 0.1, u_mid + 0.9);
    let b = quad_face(
        &mut brep,
        [
            surface_b.point_at(u0, x_range.0),
            surface_b.point_at(u0, x_range.1),
            surface_b.point_at(u1, x_range.1),
            surface_b.point_at(u1, x_range.0),
        ],
        surface_b,
    );
    (brep, a, b)
}

#[test]
fn crossing_cylinder_walls_are_flagged() {
    let (brep, a, b) = two_walls((0.1, 0.6));
    assert!(
        matches!(
            brep.validate_self_intersection_among(&[a, b]),
            Err(ValidationError::SelfIntersection { .. })
        ),
        "the walls cross along x^2 + y^2 = 0.25"
    );
}

#[test]
fn cylinder_walls_apart_are_not_flagged() {
    let (brep, a, b) = two_walls((0.6, 1.0));
    let verdict = brep.validate_self_intersection_among(&[a, b]);
    assert!(verdict.is_ok(), "{verdict:?}");
}

#[test]
fn a_near_sample_refines_onto_both_cylinders_or_is_a_near_miss() {
    let a = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.5,
    };
    let b = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.5),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 0.4,
    };
    let tol = 1e-6;
    let crossing = super::self_intersection_pair_fast_paths::cylinder_pair_crossing_near(
        &a,
        &b,
        Point3::new(0.31, 0.39, 0.52),
        tol,
    )
    .expect("the sheets cross nearby");
    assert!(((crossing.x.powi(2) + crossing.y.powi(2)).sqrt() - 0.5).abs() <= tol);
    let (dy, dz) = (crossing.y, crossing.z - 0.5);
    assert!(((dy * dy + dz * dz).sqrt() - 0.4).abs() <= tol);

    // Parallel sheets 5 tolerances apart never meet.
    let near = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.5 + 5.0 * tol,
    };
    assert!(
        super::self_intersection_pair_fast_paths::cylinder_pair_crossing_near(
            &a,
            &near,
            Point3::new(0.5, 0.0, 0.3),
            tol,
        )
        .is_none()
    );
}
