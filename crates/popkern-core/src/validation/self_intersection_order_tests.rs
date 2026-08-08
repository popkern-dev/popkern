//! PK-251 regressions: pair-order insensitivity of the
//! self-intersection verdict and the interior-overlap samples API.

use crate::error::ValidationError;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Degree-1 bilinear NURBS quad (duplicated from the sibling test
/// module; both stay under the file-size budget).
fn bilinear_patch(p00: Point3, p01: Point3, p10: Point3, p11: Point3) -> SurfaceGeometry {
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p00, p01, p10, p11],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    })
}

/// PK-251 regression: the pair verdict must not depend on which face
/// comes first in `face_keys`. The intersection solver is
/// order-sensitive for free-form×analytic pairs (measured on the
/// mixed-radius fillet corners: cylinder-first found the strip×cap
/// crossing, NURBS-first returned nothing), so the check canonicalizes
/// each pair analytic-first before intersecting.
fn crossing_cylinder_and_nurbs_fixture() -> (
    BRep,
    crate::topology::FaceHandle,
    crate::topology::FaceHandle,
) {
    let mut brep = BRep::new(Tolerance::default());

    // Face A — quarter cylinder wall, radius 0.5, axis +Z, z in [0,1].
    let mvfs_a = brep.mvfs(Point3::new(0.5, 0.0, 0.0)).unwrap();
    let mev_a1 = brep
        .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(0.0, 0.5, 0.0))
        .unwrap();
    let mev_a2 = brep
        .mev(mev_a1.half_edge_twin, Point3::new(0.0, 0.5, 1.0))
        .unwrap();
    let mev_a3 = brep
        .mev(mev_a2.half_edge_twin, Point3::new(0.5, 0.0, 1.0))
        .unwrap();
    let mef_a = brep.mef(mev_a3.half_edge_twin, mev_a1.half_edge).unwrap();
    brep.faces[mef_a.face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 0.5,
        },
        orientation: true,
    });

    // Face B — flat horizontal quad at z=0.5 crossing the wall.
    let mvfs_b = brep.mvfs(Point3::new(0.1, 0.1, 0.5)).unwrap();
    let mev_b1 = brep
        .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.1, 0.5))
        .unwrap();
    let mev_b2 = brep
        .mev(mev_b1.half_edge_twin, Point3::new(1.0, 1.0, 0.5))
        .unwrap();
    let mev_b3 = brep
        .mev(mev_b2.half_edge_twin, Point3::new(0.1, 1.0, 0.5))
        .unwrap();
    let mef_b = brep.mef(mev_b3.half_edge_twin, mev_b1.half_edge).unwrap();
    brep.faces[mef_b.face].surface = Some(SurfaceOnFace {
        geometry: bilinear_patch(
            Point3::new(0.1, 0.1, 0.5),
            Point3::new(0.1, 1.0, 0.5),
            Point3::new(1.0, 0.1, 0.5),
            Point3::new(1.0, 1.0, 0.5),
        ),
        orientation: true,
    });

    (brep, mef_a.face, mef_b.face)
}

#[test]
fn nurbs_analytic_pair_verdict_is_order_insensitive() {
    let (brep, cylinder_face, nurbs_face) = crossing_cylinder_and_nurbs_fixture();
    let analytic_first = brep.validate_self_intersection_among(&[cylinder_face, nurbs_face]);
    let free_form_first = brep.validate_self_intersection_among(&[nurbs_face, cylinder_face]);
    assert!(
        matches!(
            analytic_first,
            Err(ValidationError::SelfIntersection { .. })
        ),
        "analytic-first order must flag the crossing: {analytic_first:?}"
    );
    assert!(
        matches!(
            free_form_first,
            Err(ValidationError::SelfIntersection { .. })
        ),
        "free-form-first order must flag the same crossing: {free_form_first:?}"
    );
}

#[test]
fn interior_overlap_samples_report_both_faces_uvs_in_caller_order() {
    let (brep, cylinder_face, nurbs_face) = crossing_cylinder_and_nurbs_fixture();
    let samples = brep.interior_overlap_samples(nurbs_face, cylinder_face);
    assert!(
        !samples.is_empty(),
        "crossing pair must produce interior overlap samples"
    );
    for sample in &samples {
        // Caller order: uv_a belongs to the NURBS quad (unit domain),
        // uv_b to the cylinder chart (angle, height).
        assert!(
            (-0.01..=1.01).contains(&sample.uv_a.0) && (-0.01..=1.01).contains(&sample.uv_a.1),
            "uv_a must live in the NURBS patch domain: {:?}",
            sample.uv_a
        );
        assert!(
            (sample.point.z - 0.5).abs() < 1.0e-6,
            "crossing samples lie in the quad plane: {:?}",
            sample.point
        );
        let radial = (sample.point.x * sample.point.x + sample.point.y * sample.point.y).sqrt();
        assert!(
            (radial - 0.5).abs() < 1.0e-3,
            "crossing samples lie on the cylinder: radial={radial}"
        );
    }
}
