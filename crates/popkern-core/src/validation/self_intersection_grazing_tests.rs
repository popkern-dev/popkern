//! PK-230 pins for the material-two-sided crossing rule in the
//! analytic × free-form fast path: tangential contact whose
//! penetration never exceeds the linear tolerance is NOT a
//! self-intersection (measured root: cube_extreme's cap edge dips
//! 5.3e-9 below its support plane at the shared setback vertex and
//! was flagged as an interior crossing), while material penetration
//! past the tolerance must keep flagging.

use crate::error::ValidationError;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Plane face over [0,1]² at z=0 plus a disconnected quadratic-arch
/// NURBS face whose (u0,v0) corner sits at (0.5, 0.5, corner_z) —
/// touching the plane's INTERIOR — and rises steeply away from it.
fn plane_and_corner_touch_patch(corner_z: f64) -> BRep {
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);

    // Face A: unit quad carrying an analytic Plane at z = 0.
    let mvfs_a = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev_a1 = brep
        .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev_a2 = brep
        .mev(mev_a1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev_a3 = brep
        .mev(mev_a2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef_a = brep.mef(mev_a3.half_edge_twin, mev_a1.half_edge).unwrap();
    brep.faces[mef_a.face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    });

    // Face B: disconnected quad whose surface is a degree-(2,1) arch
    // touching the plane interior at its (u0,v0) corner and rising to
    // z = 0.4 elsewhere.
    let c = corner_z;
    let mvfs_b = brep.mvfs(Point3::new(0.5, 0.5, c)).unwrap();
    let mev_b1 = brep
        .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.5, 0.3))
        .unwrap();
    let mev_b2 = brep
        .mev(mev_b1.half_edge_twin, Point3::new(1.0, 1.0, 0.4))
        .unwrap();
    let mev_b3 = brep
        .mev(mev_b2.half_edge_twin, Point3::new(0.5, 1.0, 0.15))
        .unwrap();
    let mef_b = brep.mef(mev_b3.half_edge_twin, mev_b1.half_edge).unwrap();
    brep.faces[mef_b.face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(NurbsSurface {
            degree_u: 2,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            // Row-major: u varies within each row per NurbsSurface's
            // convention (num_u columns × num_v rows).
            control_points: vec![
                Point3::new(0.5, 0.5, c),
                Point3::new(0.75, 0.5, 0.1),
                Point3::new(1.0, 0.5, 0.3),
                Point3::new(0.5, 1.0, 0.15),
                Point3::new(0.75, 1.0, 0.25),
                Point3::new(1.0, 1.0, 0.4),
            ],
            num_u: 3,
            num_v: 2,
            weights: vec![1.0; 6],
        }),
        orientation: true,
    });

    brep
}

/// Sub-tolerance corner graze on the plane's interior: the patch
/// corner dips 1e-9 below the plane (three orders under the 1e-6
/// linear tolerance — the cube_extreme signature). The exact-zero
/// `sd <= 0` sign flip used to bisect this into samples that
/// classified (plane Interior, patch Boundary) and flagged; the
/// material-two-sided rule clears it.
#[test]
fn sub_tolerance_corner_graze_is_not_a_crossing() {
    let brep = plane_and_corner_touch_patch(-1.0e-9);
    assert!(
        brep.validate_self_intersection().is_ok(),
        "a 1e-9 corner graze must not be flagged as self-intersection"
    );
}

/// Material counterpart: the same patch sunk 2e-5 below the plane
/// (20× the linear tolerance) crosses it inside both trims and must
/// keep flagging — the honesty rule must not over-silence.
#[test]
fn material_corner_penetration_still_flags() {
    let brep = plane_and_corner_touch_patch(-2.0e-5);
    assert!(
        matches!(
            brep.validate_self_intersection(),
            Err(ValidationError::SelfIntersection { .. })
        ),
        "a 2e-5 corner penetration must keep flagging"
    );
}
