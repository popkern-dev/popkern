//! PK-669 pin for classification against a chamfered hex body: an
//! M6 nut workpiece (hex blank, both 30° corner chamfers, plain bore)
//! must classify bore-void points Outside and wall points Inside,
//! although some rays leave through a flat just beside a chamfer
//! edge, where the chamfer cone's extension lies a hair outside the
//! hex (the phantom cone hit counted a second crossing).

use popkern_bool::{PointClassification, point_in_solid_final};
use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;

use crate::builders::{hex_corner_ring, hexagon};
use crate::engine::{extrude_shape, subtract};

#[test]
fn chamfered_nut_workpiece_classifies_void_and_wall_points() {
    let tol = Tolerance::default();
    let (s, m, r_bore) = (10.0, 5.2, 2.6088);
    let mut body = extrude_shape(hexagon(s), 0.0, m, "blank", tol).expect("blank");
    for (z, outward) in [(0.0, -1.0), (m, 1.0)] {
        let ring = hex_corner_ring(s, z, outward, tol).expect("ring");
        body = subtract(&body, &ring, "chamfer", tol).expect("chamfer");
    }
    let bore = extrude_shape(
        popkern_sketch::ProfileShape::Circle { radius: r_bore },
        -0.1 * m,
        1.2 * m,
        "bore",
        tol,
    )
    .expect("bore");
    body = subtract(&body, &bore, "bore", tol).expect("bored");
    let aabbs: Vec<_> = body
        .faces
        .keys()
        .filter_map(|fh| {
            let pts = body.face_boundary_sample_points(fh, 32);
            Aabb::from_points(&pts).map(|a| (fh, a))
        })
        .collect();
    let mut wrong = Vec::new();
    for (r, expected) in [
        (0.0, PointClassification::Outside),
        (2.0, PointClassification::Outside),
        (4.0, PointClassification::Inside),
    ] {
        for z in [0.18, 0.3, 2.6, 4.9, 5.02] {
            for k in 0..16 {
                let a = std::f64::consts::TAU * k as f64 / 16.0;
                let p = Point3::new(r * a.cos(), r * a.sin(), z);
                let got = point_in_solid_final(&p, &body, &aabbs, &tol).expect("classifies");
                if got != expected {
                    wrong.push(format!("r={r} z={z} k={k}: {got:?}"));
                }
            }
        }
    }
    assert!(wrong.is_empty(), "misclassified: {wrong:?}");
}
