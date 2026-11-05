// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::*;

fn square(half: f64, z: f64) -> Vec<Point3> {
    vec![
        Point3::new(-half, -half, z),
        Point3::new(half, -half, z),
        Point3::new(half, half, z),
        Point3::new(-half, half, z),
    ]
}

/// The enclosed volume: Σ area · (outward normal · corner) / 3 over the
/// planar faces, the outward normal the surface's own (with the face's
/// orientation). Positive, and the solid's volume, only when every face
/// points out of the material.
fn signed_volume(brep: &BRep) -> f64 {
    brep.faces
        .values()
        .zip(brep.faces.keys())
        .map(|(face, fh)| {
            let Some(sof) = face.surface.as_ref() else {
                return f64::NAN;
            };
            let crate::geometry::surface::SurfaceGeometry::Plane { normal, .. } = sof.geometry
            else {
                return f64::NAN;
            };
            let outward = if sof.orientation { normal } else { -normal };
            let loop_ = brep.face_outer_loop_positions(fh);
            let area = newell(&loop_).length() * 0.5;
            area * outward.dot(&(loop_[0] - Point3::new(0.0, 0.0, 0.0))) / 3.0
        })
        .sum()
}

/// The Euler prism is built inward (−4); turned outward it measures
/// +4 and validates, and turning an outward body is a no-op.
#[test]
fn an_inward_prism_is_turned_outward() {
    let mut prism = BRep::prism(
        &square(1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .expect("prism");
    assert!((signed_volume(&prism) + 4.0).abs() < 1e-12);
    assert!(prism.orient_planar_outward().expect("planar"));
    assert!((signed_volume(&prism) - 4.0).abs() < 1e-12);
    prism.validate().expect("turned prism validates");
    assert!(!prism.orient_planar_outward().expect("planar"));
}

#[test]
fn a_square_frustum_is_a_prism_topology_facing_out() {
    let brep = BRep::polygon_frustum(&square(2.0, 0.0), &square(1.0, 3.0), Tolerance::default())
        .expect("frustum");
    assert_eq!(
        (brep.vertices.len(), brep.edges.len(), brep.faces.len()),
        (8, 12, 6)
    );
    // h/3 (A1 + A2 + sqrt(A1 A2)) = 1 (16 + 4 + 8).
    let v = signed_volume(&brep);
    assert!((v - 28.0).abs() < 1e-12, "{v}");
}

#[test]
fn either_winding_faces_out() {
    let bottom: Vec<Point3> = square(2.0, 0.0).into_iter().rev().collect();
    let top: Vec<Point3> = square(1.0, 3.0).into_iter().rev().collect();
    let brep = BRep::polygon_frustum(&bottom, &top, Tolerance::default()).expect("frustum");
    let v = signed_volume(&brep);
    assert!((v - 28.0).abs() < 1e-12, "{v}");
}

/// A top turned against the bottom twists every side: declined, not
/// served with fitted planes.
#[test]
fn a_twisted_side_declines() {
    let mut top = square(1.0, 3.0);
    top.rotate_left(1);
    let err =
        BRep::polygon_frustum(&square(2.0, 0.0), &top, Tolerance::default()).expect_err("twisted");
    assert!(err.to_string().contains("not planar"), "{err}");
}

#[test]
fn a_flat_or_mismatched_pair_declines() {
    let flat = BRep::polygon_frustum(&square(2.0, 0.0), &square(1.0, 0.0), Tolerance::default())
        .expect_err("flat");
    assert!(flat.to_string().contains("stand off"), "{flat}");
    let short = BRep::polygon_frustum(
        &square(2.0, 0.0),
        &square(1.0, 3.0)[..3],
        Tolerance::default(),
    )
    .expect_err("mismatched");
    assert!(short.to_string().contains("same corner count"), "{short}");
}
