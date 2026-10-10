//! Edge-face detection on reversed arcs, in-surface edges and endpoint
//! touches (PK-750).

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{AnalyticProfile, ProfileShape, SketchPlane};

use crate::intersect::compute_face_aabbs;
use crate::intersect_edge_face::detect_edge_face_intersections;

/// A sketch rod (radius 4.5, height 16) on the z axis whose circle
/// starts at +y: its rims are two half-circle edges each, the top ones
/// stored with t_start > t_end.
fn sketch_rod() -> BRep {
    let sheet = AnalyticProfile {
        shape: ProfileShape::Circle { radius: 4.5 },
        center: [0.0, 0.0],
    }
    .to_brep(
        SketchPlane {
            origin: Point3::new(0.0, 0.0, 0.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(-1.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        },
        Tolerance::default(),
        None,
    )
    .expect("sheet");
    let face = sheet.regions[0].face;
    let mut b = sheet.brep;
    b.extrude(face, Vec3::new(0.0, 0.0, 1.0), 16.0, None)
        .expect("extrude");
    b
}

#[test]
fn reversed_rim_arcs_are_crossed_where_a_slab_cuts_them() {
    let rod = sketch_rod();
    let reversed = rod
        .edges
        .values()
        .filter_map(|e| e.curve.as_ref())
        .filter(|c| c.t_start > c.t_end)
        .count();
    assert!(reversed > 0, "the fixture must carry reversed arcs");
    // A slab x ∈ [−1, 1] across the top rim (z ∈ [10, 20]).
    let profile = [
        Point3::new(-1.0, -10.0, 10.0),
        Point3::new(1.0, -10.0, 10.0),
        Point3::new(1.0, 10.0, 10.0),
        Point3::new(-1.0, 10.0, 10.0),
    ];
    let slab = BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        10.0,
        Tolerance::default(),
    )
    .expect("slab");
    let tol = Tolerance::default();
    let hits = detect_edge_face_intersections(
        &rod,
        &slab,
        &compute_face_aabbs(&slab),
        &compute_face_aabbs(&rod),
        &tol,
    );
    // The top rim crosses x = ±1 at y = ±√(4.5² − 1): four points.
    let rim: Vec<_> = hits
        .iter()
        .filter(|h| h.edge_of_a && (h.point.z - 16.0).abs() < 1e-9)
        .collect();
    assert_eq!(rim.len(), 4, "{rim:#?}");
    for h in rim {
        let r = h.point.x.hypot(h.point.y);
        assert!(
            (r - 4.5).abs() < 1e-9 && (h.point.x.abs() - 1.0).abs() < 1e-9,
            "{h:?}"
        );
    }
}
