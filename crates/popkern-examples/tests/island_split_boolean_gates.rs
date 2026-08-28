//! PK-341: island-split fragments on faces with concentric rings.
//!
//! `cylinder − tube` island-splits the tube's annular caps along a
//! concentric interior circle. Before the fix the 3D-Newell walker
//! election called both circle copies "separate", hole nesting
//! rejected the concentric inner rim on winding, classification
//! sampled knife-edge points ON the other operand, and the served
//! body was NON-CLOSED yet validate_all-clean. Now the chart-role
//! path owns island-pair faces (σ-folded roles + provenance hole
//! override), classification samples step into fragment interiors,
//! and the ray cast clamps plane hits to boundary boxes — the
//! subtractive construction is topology-identical to the additive
//! stacked-union construction.

use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::grooved_cylinder;
use popkern_mesh::{TessellationOptions, tessellate};
use std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn mesh_stats(brep: &BRep) -> f64 {
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
fn subtractive_grooved_cylinder_matches_the_additive_construction() {
    let base = BRep::cylinder(0.5, 1.0, tol()).unwrap();
    let ring = BRep::cylinder(0.7, 0.1, tol())
        .unwrap()
        .transform(&Transform3::translation(0.0, 0.0, 0.45))
        .unwrap();
    let plug = BRep::cylinder(0.3, 0.3, tol())
        .unwrap()
        .transform(&Transform3::translation(0.0, 0.0, 0.35))
        .unwrap();
    let tube = popkern_bool::subtraction(ring, plug, tol(), None).expect("tube");
    assert!(tube.validate_all().is_empty());
    let result = popkern_bool::subtraction(base, tube, tol(), None).expect("grooved");

    // Exact grooved topology — identical to the additive fixture:
    // 2 wall bands + groove floor (cylinders), 2 caps + 2 annular
    // groove walls WITH their inner rims (planes).
    assert_eq!(result.faces.len(), 7, "faces");
    assert_eq!(result.edges.len(), 9, "edges");
    assert_eq!(result.vertices.len(), 6, "vertices");
    let holed: usize = result
        .faces
        .values()
        .filter(|f| f.inner_loops.len() == 1)
        .count();
    assert_eq!(holed, 2, "both groove walls keep their inner rims");
    assert!(result.validate_all().is_empty());

    // Grooved closed-form volume (position-welded meshes make the
    // signed volume the sealed-ness witness).
    let volume = mesh_stats(&result);
    let expected = PI * (0.25 - 0.016);
    let rel = (volume - expected).abs() / expected;
    assert!(
        rel < 0.01,
        "volume {volume:.5} vs {expected:.5} (rel {rel:.5})"
    );

    // And it matches the additive corpus fixture's census exactly.
    let additive = grooved_cylinder(tol()).unwrap();
    assert_eq!(result.faces.len(), additive.expected_faces);
    assert_eq!(result.edges.len(), additive.expected_edges);
    assert_eq!(result.vertices.len(), additive.expected_vertices);
    assert_eq!(
        popkern_examples::offset_corpus::surface_census(&result),
        additive.expected_census
    );
}
