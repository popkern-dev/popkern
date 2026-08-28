//! PK-345 gate: `point_in_solid_final` classifies probe points nudged
//! across every corner sphere of the filleted cube correctly — the
//! padded parameter box covers the whole sphere, and the phantom 7/8
//! flipped 4 of the 8 identical corner probes before the sphere-cap
//! hit gate.

use popkern_bool::{PointClassification, point_in_solid_final};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_examples::offset_corpus::filleted_cube;

#[test]
fn corner_probes_classify_correctly_on_the_filleted_cube() {
    let f = filleted_cube(Tolerance::default()).expect("fixture");
    // The classifier's broad-phase AABBs are internal; rebuild the
    // equivalent sampled boxes from the public surface data.
    let aabbs: Vec<_> = f
        .brep
        .faces
        .keys()
        .filter_map(|fh| {
            let positions = f.brep.face_boundary_sample_points(fh, 8);
            let sof = f.brep.faces[fh].surface.as_ref()?;
            use popkern_core::geometry::surface::SurfaceEval;
            use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
            let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&positions, &sof.geometry);
            let mut pts = Vec::new();
            for i in 0..=4 {
                for j in 0..=4 {
                    let u = u0 + (u1 - u0) * i as f64 / 4.0;
                    let v = v0 + (v1 - v0) * j as f64 / 4.0;
                    pts.push(sof.geometry.point_at(u, v));
                }
            }
            popkern_core::geometry::intersection::surface_subdivision::Aabb::from_points(&pts)
                .map(|a| (fh, a))
        })
        .collect();
    let tol = Tolerance::default();
    let eps = 2.0e-4;
    for cx in [0.1, 0.9] {
        for cy in [0.1, 0.9] {
            for cz in [0.1, 0.9] {
                let c = Point3::new(cx, cy, cz);
                let dir = Vec3::new(
                    if cx < 0.5 { -1.0 } else { 1.0 },
                    if cy < 0.5 { -1.0 } else { 1.0 },
                    if cz < 0.5 { -1.0 } else { 1.0 },
                ) * (1.0 / 3f64.sqrt());
                let p = c + dir * 0.1;
                let outside = point_in_solid_final(&(p + dir * eps), &f.brep, &aabbs, &tol)
                    .expect("probe classifies");
                let inside = point_in_solid_final(&(p - dir * eps), &f.brep, &aabbs, &tol)
                    .expect("probe classifies");
                assert_eq!(
                    outside,
                    PointClassification::Outside,
                    "corner ({cx},{cy},{cz}): outward probe"
                );
                assert_eq!(
                    inside,
                    PointClassification::Inside,
                    "corner ({cx},{cy},{cz}): inward probe"
                );
            }
        }
    }
}
