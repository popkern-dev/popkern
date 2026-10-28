//! Reflections through `BRep::transform` (PK-976).

use crate::geometry::surface::SurfaceEval;
use crate::geometry::transform::Transform3;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// A reflection keeps rotational faces honest: their derived chart
/// basis stays right-handed, so a mirrored body's pcurves must reflect
/// (u -> delta - u) and their orientation flags toggle, or a cone
/// mirrored for a nut's lower chamfer keeps shifted pcurves and a flag
/// folding its normal outward beside its inward-folded base.
#[test]
fn reflection_keeps_rotational_faces_honest() {
    let tol = Tolerance::default();
    let mirror = Transform3::reflection(
        Point3::new(0.3, -0.2, 1.1),
        Vec3::new(0.2, -0.4, 1.0).normalized().unwrap(),
    )
    .unwrap();
    // (name, body, a point inside it: every folded normal of a convex
    // body points to one side of it).
    let bodies = [
        (
            "cylinder",
            BRep::cylinder(1.5, 2.0, tol).unwrap(),
            Some(Point3::new(0.0, 0.0, 1.0)),
        ),
        (
            "cone",
            BRep::cone(2.0, 3.0, tol).unwrap(),
            Some(Point3::new(0.0, 0.0, 0.5)),
        ),
        (
            "sphere",
            BRep::sphere(1.25, tol).unwrap(),
            Some(Point3::new(0.0, 0.0, 0.0)),
        ),
        ("torus", BRep::torus(3.0, 1.0, tol).unwrap(), None),
    ];
    for (name, body, inside) in bodies {
        let mirrored = body.transform(&mirror).unwrap();
        assert!(
            mirrored.dishonest_winding_faces().is_empty(),
            "{name}: mirrored faces wind against their flags"
        );
        for (heh, _) in &mirrored.half_edges {
            if let Some(d) = mirrored.measure_half_edge_pcurve_deviation(heh) {
                assert!(
                    d.max_deviation < 1e-6,
                    "{name}: pcurve of {heh:?} sits {} off its mirrored edge",
                    d.max_deviation
                );
            }
        }
        let Some(inside) = inside else {
            continue;
        };
        let centre = mirror.apply_point(&inside);
        let mut sides = Vec::new();
        for (fh, face) in &mirrored.faces {
            let sof = face.surface.as_ref().unwrap();
            let samples = mirrored.face_boundary_sample_points(fh, 8);
            // Off the boundary: towards the centre from a boundary
            // sample, pulled back onto the surface.
            let a = samples[0];
            let toward = Point3::new(
                a.x + 0.1 * (centre.x - a.x),
                a.y + 0.1 * (centre.y - a.y),
                a.z + 0.1 * (centre.z - a.z),
            );
            let (u, v) = sof.geometry.closest_point(&toward).unwrap();
            let Some(n) = sof.folded_normal(u, v) else {
                continue;
            };
            sides.push(n.dot(&(sof.geometry.point_at(u, v) - centre)) > 0.0);
        }
        assert!(
            sides.windows(2).all(|w| w[0] == w[1]),
            "{name}: mirrored faces fold their normals to different sides: {sides:?}"
        );
    }
}
