//! PK-707: ray classification against cylinders longer than the
//! historical ±100 parameter window. Every lateral hit beyond |v| = 100
//! was filtered, so points inside a long bore classified Outside (the
//! cellular boolean then kept a 300-long piano hinge leaf's bore disc
//! and fell back; popkern-hinge's tests/piano.rs pins that case).

use popkern_bool::{PointClassification, point_in_solid_final};
use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn classify(p: Point3, body: &BRep) -> PointClassification {
    let aabbs: Vec<_> = body
        .faces
        .keys()
        .filter_map(|fh| {
            Aabb::from_points(&body.face_boundary_sample_points(fh, 32)).map(|a| (fh, a))
        })
        .collect();
    point_in_solid_final(&p, body, &aabbs, &tol()).expect("classifies")
}

#[test]
fn points_inside_a_long_cylinder_classify_inside() {
    let rod = BRep::cylinder(1.4, 300.0, tol()).expect("rod");
    for z in [5.0, 99.0, 150.0, 250.0, 299.0] {
        for (x, y) in [(0.0, -0.885), (0.3, 0.2), (-1.0, 0.5)] {
            assert_eq!(
                classify(Point3::new(x, y, z), &rod),
                PointClassification::Inside,
                "({x}, {y}, {z})"
            );
        }
        assert_eq!(
            classify(Point3::new(2.0, 0.0, z), &rod),
            PointClassification::Outside
        );
    }
}
