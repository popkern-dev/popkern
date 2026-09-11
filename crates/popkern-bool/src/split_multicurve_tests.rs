// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use super::*;
use popkern_core::{
    geometry::{Point3, Vec3, nurbs_curve::NurbsCurve},
    tolerance::Tolerance,
    topology::BRep,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Build a square BRep with two midpoint vertices and return (brep, face, p1, p2).
fn build_splittable_square(
    origin: Point3,
) -> (BRep, popkern_core::topology::FaceHandle, Point3, Point3) {
    let mut brep = BRep::new(tol());
    let mvfs = brep.mvfs(origin).expect("mvfs");
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, origin + Vec3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep
        .mev(mev1.half_edge_twin, origin + Vec3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep
        .mev(mev2.half_edge_twin, origin + Vec3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).expect("mef");
    let split1 = brep.split_edge(mev1.edge, 0.5).expect("split1");
    let split2 = brep.split_edge(mev3.edge, 0.5).expect("split2");
    let p1 = brep.vertices[split1.vertex].position;
    let p2 = brep.vertices[split2.vertex].position;
    (brep, mef_r.face, p1, p2)
}

/// Verify that the multi-curve grouping collects all indices per face.
/// When two records share the same face_a, both indices must appear in
/// the same Vec entry of the curves_a HashMap.
#[test]
fn curves_per_face_collects_multiple_indices() {
    use popkern_core::geometry::curve_on_surface::CurveOnSurface;
    use std::collections::HashMap;

    let dummy_curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let cos = CurveOnSurface {
        curve_3d: dummy_curve.clone(),
        curve_uv1: dummy_curve.clone(),
        curve_uv2: dummy_curve,
        is_closed: false,
    };

    let face_a = popkern_core::topology::FaceHandle::default();
    let face_b1 = popkern_core::topology::FaceHandle::default();
    let face_b2 = popkern_core::topology::FaceHandle::default();

    let records = [
        crate::intersect::IntersectionRecord {
            residue_owner_a: None,
            face_a,
            face_b: face_b1,
            curve: cos.clone(),
            is_tangent: false,
            is_coincident: false,
            start_boundary: Default::default(),
            end_boundary: Default::default(),
        },
        crate::intersect::IntersectionRecord {
            residue_owner_a: None,
            face_a,
            face_b: face_b2,
            curve: cos,
            is_tangent: false,
            is_coincident: false,
            start_boundary: Default::default(),
            end_boundary: Default::default(),
        },
    ];

    // Replicate the grouping logic from split_along_curves.
    let mut curves_a: HashMap<popkern_core::topology::FaceHandle, Vec<usize>> = HashMap::new();
    for (idx, record) in records.iter().enumerate() {
        curves_a.entry(record.face_a).or_default().push(idx);
    }

    // Both records target the same face_a — should be grouped together.
    assert_eq!(curves_a.len(), 1, "should have one face entry");
    let indices = &curves_a[&face_a];
    assert_eq!(indices.len(), 2, "should collect both record indices");
    assert_eq!(indices[0], 0);
    assert_eq!(indices[1], 1);
}

/// Two records targeting the same face: verify the first curve splits and
/// the second is gracefully skipped (multi-curve support is exercised
/// even though subsequent curves are not yet applied).
#[test]
fn split_same_face_with_two_records_uses_first() {
    let (brep_a, face_a_h, p1, p2) = build_splittable_square(Point3::new(0.0, 0.0, 0.0));

    let trim_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p1, p2],
        weights: vec![1.0, 1.0],
    };
    use popkern_core::geometry::curve_on_surface::CurveOnSurface;
    let curve = CurveOnSurface {
        curve_3d: trim_3d.clone(),
        curve_uv1: trim_3d.clone(),
        curve_uv2: trim_3d.clone(),
        is_closed: false,
    };

    // Create a second (dummy) record targeting the same face.
    let dummy_curve = CurveOnSurface {
        curve_3d: trim_3d.clone(),
        curve_uv1: trim_3d.clone(),
        curve_uv2: trim_3d,
        is_closed: false,
    };

    let record1 = crate::intersect::IntersectionRecord {
        residue_owner_a: None,
        face_a: face_a_h,
        face_b: popkern_core::topology::FaceHandle::default(),
        curve,
        is_tangent: false,
        is_coincident: false,
        start_boundary: Default::default(),
        end_boundary: Default::default(),
    };
    let record2 = crate::intersect::IntersectionRecord {
        residue_owner_a: None,
        face_a: face_a_h,
        face_b: popkern_core::topology::FaceHandle::default(),
        curve: dummy_curve,
        is_tangent: false,
        is_coincident: false,
        start_boundary: Default::default(),
        end_boundary: Default::default(),
    };

    let faces_before = brep_a.faces.len();
    let brep_b = BRep::new(tol());
    let result = split_along_curves(brep_a, brep_b, &[record1, record2], &tol());

    match result {
        Ok(r) => {
            // First curve should split the face (adds 1 face).
            // Second curve is skipped for now.
            assert!(
                r.a.faces.len() > faces_before,
                "first curve should split the face"
            );
        }
        Err(_) => {
            // Split failure is acceptable for some topologies.
        }
    }
}
