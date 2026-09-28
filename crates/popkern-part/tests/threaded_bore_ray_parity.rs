//! PK-566 gate: point-in-solid classification inside and around an
//! M6 nut's threaded bore.
//!
//! Rays cast from inside the bore cross the CREST BAND, a helical
//! strip on the bore cylinder between the thread grooves. The trim
//! test that accepts those hits ran a horizontal even-odd on the
//! band's FOLDED trim polygon. That polygon's turn pieces restart at
//! the seam without any seam segment, so wherever the strip wrapped
//! the seam at the hit's height the test read the complement. A
//! phantom crest hit then flipped the parity of three of five rays,
//! and bore-void points classified Inside. That sent the mouth
//! chamfer's void-side tool wall into the cut, and the cellular
//! assembly went non-manifold.

use popkern_bool::{PointClassification, point_in_solid_final};
use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_part::fastener::{HexNut, ThreadSpec};

/// The classifier's broad-phase AABBs are internal; rebuild
/// equivalent boxes from the public surface data, sampled densely
/// enough (16 × 512) to cover every turn of a multi-turn helical
/// flank. A coarse grid clips a flank's box and rejects genuine hits.
fn sampled_face_aabbs(brep: &BRep) -> Vec<(FaceHandle, Aabb)> {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
    brep.faces
        .keys()
        .filter_map(|fh| {
            let positions = brep.face_boundary_sample_points(fh, 8);
            let sof = brep.faces[fh].surface.as_ref()?;
            let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&positions, &sof.geometry);
            let mut pts = Vec::new();
            for i in 0..=16 {
                for j in 0..=512 {
                    let u = u0 + (u1 - u0) * i as f64 / 16.0;
                    let v = v0 + (v1 - v0) * j as f64 / 512.0;
                    pts.push(sof.geometry.point_at(u, v));
                }
            }
            Aabb::from_points(&pts).map(|a| (fh, a))
        })
        .collect()
}

#[test]
fn threaded_bore_void_and_wall_points_classify_correctly() {
    let tol = Tolerance::default();
    let blank = HexNut::new(ThreadSpec::parse("M6").expect("M6"))
        .clearance(0.15)
        .mouth_chamfer(false)
        .build(tol)
        .expect("threaded blank");
    let nut = blank.brep();
    let aabbs = sampled_face_aabbs(nut);
    let mut wrong = Vec::new();
    // r 2.0 and 2.4837 lie inside the 2.6087 bore (void). r 4.0 lies
    // between the ~3.2 groove base and the 5.0 hex flats (solid).
    for (r, expected) in [
        (2.0, PointClassification::Outside),
        (2.4837, PointClassification::Outside),
        (4.0, PointClassification::Inside),
    ] {
        for z in [0.18, 1.0, 2.6] {
            for k in 0..16 {
                let a = std::f64::consts::TAU * k as f64 / 16.0;
                let p = Point3::new(r * a.cos(), r * a.sin(), z);
                let got = point_in_solid_final(&p, nut, &aabbs, &tol).expect("classifies");
                if got != expected {
                    wrong.push(format!("r={r} z={z} k={k}: {got:?}"));
                }
            }
        }
    }
    assert!(wrong.is_empty(), "misclassified: {wrong:?}");
}
