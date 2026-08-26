//! PK-343: cone-apex loop closure through booleans.
//!
//! The legacy walk refused the apex TURNAROUND (a cone lateral's seam
//! halves are twins inside the same loop; U-turn prevention barred the
//! only continuation), dropping the seam-to-apex piece; the election's
//! size rule then filed the small apex region as a "hole" of the
//! disjoint bottom frustum. Every cone-through-body union served a
//! torn or apex-less result until the PK-342 closure gate exposed it.

use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{TessellationOptions, tessellate};
use std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn volume(brep: &BRep) -> f64 {
    let mesh = tessellate(brep, &TessellationOptions::default()).expect("tessellate");
    let mut vol = 0.0_f64;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (vol / 6.0).abs()
}

#[test]
fn cone_through_cube_union_serves_the_apex_whole() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let cone = BRep::cone(0.4, 2.0, tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let result = popkern_bool::union(cube, cone, tol(), None).expect("union serves");

    // 6 cube faces (top/bottom holed) + frustum lateral + base cap +
    // apex piece.
    assert_eq!(result.faces.len(), 9, "faces");
    let holed: usize = result
        .faces
        .values()
        .filter(|f| f.inner_loops.len() == 1)
        .count();
    assert!(holed >= 2, "cube top and bottom carry the cone holes");
    assert!(result.validate_all().is_empty());

    // Exact closed form: cube + below-frustum (r 0.4→0.3 over 0.5) +
    // apex tip (r 0.1, h 0.5). Cone r(z) = 0.4·(1.5 − z)/2.
    let frustum = PI * 0.5 / 3.0 * (0.4 * 0.4 + 0.4 * 0.3 + 0.3 * 0.3);
    let tip = PI * 0.1 * 0.1 * 0.5 / 3.0;
    let expected = 1.0 + frustum + tip;
    let v = volume(&result);
    let rel = (v - expected).abs() / expected;
    assert!(rel < 0.01, "volume {v:.5} vs {expected:.5} (rel {rel:.5})");
}
