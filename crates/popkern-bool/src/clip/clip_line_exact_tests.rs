use super::*;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

/// A planar cap of a radius-`r` cylinder and its height.
fn cap(r: f64) -> (BRep, FaceHandle, f64) {
    let brep = BRep::cylinder(r, 2.0, Tolerance::default()).expect("cylinder");
    let (face, z) = brep
        .faces
        .iter()
        .find_map(|(h, f)| match &f.surface.as_ref()?.geometry {
            SurfaceGeometry::Plane { origin, .. } => {
                let z = f.surface.as_ref()?.geometry.point_at(0.0, 0.0).z;
                let _ = origin;
                Some((h, z))
            }
            _ => None,
        })
        .expect("a planar cap");
    (brep, face, z)
}

#[test]
fn a_chord_of_a_circular_cap_ends_on_the_circle() {
    // The face has one boundary vertex: its vertex polygon cannot clip.
    for (r, y) in [(1.0, 0.3), (54.0, 47.7), (35.25, -32.0), (1.0, -0.999)] {
        let (brep, face, z) = cap(r);
        let c = line(Point3::new(-3.0 * r, y, z), Point3::new(3.0 * r, y, z));
        let (s0, s1) = clip_line_to_planar_face_exact(&brep, face, &c, 0.5).expect("clip");
        let half = (r * r - y * y).sqrt();
        for (s, x) in [(s0, -half), (s1, half)] {
            let p = c.point_at(s);
            assert!(
                (p.x - x).abs() < 1e-9 * r && (p.x.hypot(p.y) - r).abs() < 1e-9 * r,
                "r {r} y {y}: end ({}, {}) not on the rim at x {x}",
                p.x,
                p.y
            );
        }
    }
}

#[test]
fn a_line_missing_the_cap_or_with_its_midpoint_outside_does_not_clip() {
    let (brep, face, z) = cap(1.0);
    let miss = line(Point3::new(-3.0, 1.5, z), Point3::new(3.0, 1.5, z));
    assert!(clip_line_to_planar_face_exact(&brep, face, &miss, 0.5).is_none());
    // The line crosses the disc but its midpoint (x = 2.5) is outside.
    let off = line(Point3::new(0.0, 0.0, z), Point3::new(5.0, 0.0, z));
    assert!(clip_line_to_planar_face_exact(&brep, face, &off, 0.5).is_none());
}

#[test]
fn a_line_through_the_seam_vertex_counts_the_vertex_once() {
    // The cap's one vertex sits on the rim; a diameter through it
    // must still give one span (an even count would read outside).
    let (brep, face, z) = cap(1.0);
    let v = brep
        .vertices
        .values()
        .map(|v| v.position)
        .find(|p| (p.z - z).abs() < 1e-9)
        .expect("a rim vertex");
    let c = line(
        Point3::new(-3.0 * v.x, -3.0 * v.y, z),
        Point3::new(3.0 * v.x, 3.0 * v.y, z),
    );
    let (s0, s1) = clip_line_to_planar_face_exact(&brep, face, &c, 0.5).expect("clip");
    let (p0, p1) = (c.point_at(s0), c.point_at(s1));
    assert!(
        (p0.distance_to(&p1) - 2.0).abs() < 1e-9,
        "span {}",
        p0.distance_to(&p1)
    );
}
