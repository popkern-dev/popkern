//! PK-313 volume-envelope pin for the grazing micro-lens torus union.
//!
//! Companion to popkern-bool's `grazing_micro_lens_serves_with_exact_
//! junctions` pin: a served union can still be the WRONG SHAPE while
//! validating clean (the PK-310 uncut-overlay lesson — only a volume
//! oracle catches that class). Two tori analytic 4.8366; this
//! fixture family's tessellation measures 4.800–4.801 across the
//! whole depth ladder (discretization deficit ~0.75%), so a ±0.5%
//! band around the measured centre separates right-shape from any
//! wrong-shape class by an order of magnitude.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn micro_lens_union_volume_in_envelope() {
    if cfg!(debug_assertions) {
        eprintln!("micro-lens volume pin: skipped in debug — seconds-scale boolean");
        return;
    }
    let tol = Tolerance::default();
    for depth in [0.0005_f64, 0.002] {
        let d = 2.0 * (1.0 + 0.35) - depth;
        let ta = BRep::torus(1.0, 0.35, tol).expect("torus a");
        let tb = BRep::torus(1.0, 0.35, tol)
            .expect("torus b")
            .transform(&Transform3::translation(d, 0.0, 0.0))
            .expect("translate");
        let b = boolean(BoolOperation::Union, ta, tb, tol, None)
            .unwrap_or_else(|e| panic!("depth {depth}: micro-lens union must serve, got {e}"));
        let mesh = popkern_mesh::tessellate(&b, &popkern_mesh::TessellationOptions::default())
            .unwrap_or_else(|e| panic!("depth {depth}: tessellation failed: {e}"));
        let mut vol6 = 0.0_f64;
        for tri in &mesh.triangles {
            let a = mesh.vertices[tri[0]].position;
            let bb = mesh.vertices[tri[1]].position;
            let c = mesh.vertices[tri[2]].position;
            vol6 += a.x * (bb.y * c.z - bb.z * c.y) - a.y * (bb.x * c.z - bb.z * c.x)
                + a.z * (bb.x * c.y - bb.y * c.x);
        }
        let vol = vol6.abs() / 6.0;
        assert!(
            (vol - 4.8005).abs() < 0.024,
            "depth {depth}: mesh volume {vol} outside the measured serving envelope"
        );
    }
}
