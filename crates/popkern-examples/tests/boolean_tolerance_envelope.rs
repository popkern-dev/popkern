//! PK-342: boolean quality across the tolerance envelope.
//!
//! Two defects made coarse-tolerance booleans serve garbage:
//! the shared-node widened gate (`tolerance × 1e3`) was an ABSOLUTE
//! distance that at linear 1e-4 reached 0.2 — exactly a bore circle's
//! clearance from the cube rim — splitting rims at nodes far off them
//! (a resonance: 1e-3's even larger gate tripped its own endpoint
//! guard, masking it); and the boolean serve gate never checked loop
//! closure, shell gaps, or curve/surface agreement, so torn results
//! logged their defects and served anyway. The gate is now capped
//! geometry-relative and the closure checks decline.

use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{TessellationOptions, tessellate};

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

/// Mesh edges used by exactly one triangle (position-welded).
fn open_edges(mesh: &popkern_mesh::TriMesh) -> usize {
    use std::collections::HashMap;
    type Bits = (u64, u64, u64);
    let key = |p: popkern_core::geometry::Point3| -> Bits {
        (p.x.to_bits(), p.y.to_bits(), p.z.to_bits())
    };
    let mut count: HashMap<(Bits, Bits), usize> = HashMap::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            let a = key(mesh.vertices[t[k]].position);
            let b = key(mesh.vertices[t[(k + 1) % 3]].position);
            let e = if a <= b { (a, b) } else { (b, a) };
            *count.entry(e).or_insert(0) += 1;
        }
    }
    count.values().filter(|&&c| c != 2).count()
}

fn nurbs_edges(brep: &BRep) -> usize {
    brep.edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(popkern_core::geometry::curve::CurveGeometry::Nurbs(_))
            )
        })
        .count()
}

#[test]
fn boolean_tolerance_envelope_serves_clean_or_declines_structured() {
    let vs = 4.0 / 3.0 * std::f64::consts::PI * 0.6_f64.powi(3);
    let cap = std::f64::consts::PI * 0.01 * (3.0 * 0.6 - 0.1) / 3.0;
    let truth_sphere = 1.0 - (vs - 6.0 * cap);
    let truth_cyl = 1.0 - std::f64::consts::PI * 0.09;

    for lin in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2] {
        let t = Tolerance::default().with_linear(lin);

        // Through-cylinder subtraction: must serve clean at EVERY
        // tolerance in the envelope (the 1e-4 resonance served a
        // hole-less body with 0.2-deviated edges).
        let cube = BRep::unit_cube(t).unwrap();
        let cyl = BRep::cylinder(0.3, 1.4, t)
            .unwrap()
            .transform(&Transform3::translation(0.5, 0.5, -0.2))
            .unwrap();
        let b = popkern_bool::subtraction(cube, cyl, t, None)
            .unwrap_or_else(|e| panic!("cyl-cube at {lin:e} must serve: {e}"));
        assert!(
            b.validate_all().is_empty(),
            "cyl-cube at {lin:e} must be validator-clean"
        );
        let v = volume(&b);
        let rel = (v - truth_cyl).abs() / truth_cyl;
        assert!(
            rel < 0.01,
            "cyl-cube at {lin:e}: vol {v:.4} vs {truth_cyl:.4}"
        );

        // Torus cut by a box above z = 0.1 (off the torus's own
        // equator edges, which lie in z = 0 — a cut through them is
        // the coincident-edge class, PK-420): the kept part carries
        // two latitude rims (r = R ± r·cos(asin(0.1/r))) and the
        // pieces of the primitive's meridian edges on its torus
        // faces. Serves clean with every rim an exact arc at every
        // tolerance (PK-418); the volume below the cut is
        // 4πR ∫_{-r}^{0.1} √(r² − z²) dz, and the chordal mesh of a
        // tube under-reads it by well under 2%.
        let torus = BRep::torus(1.0, 0.3, t).unwrap();
        let box_ = BRep::unit_cube(t)
            .unwrap()
            .transform(&Transform3::scale_uniform(3.0))
            .unwrap()
            .transform(&Transform3::translation(-1.5, -1.5, 0.1))
            .unwrap();
        let truth_torus = {
            let (big_r, r, h) = (1.0_f64, 0.3_f64, 0.1_f64);
            let f = |z: f64| 0.5 * z * (r * r - z * z).sqrt() + 0.5 * r * r * (z / r).asin();
            4.0 * std::f64::consts::PI * big_r * (f(h) - f(-r))
        };
        // Serves clean at every tolerance. At 1e-2 (3% of the tube
        // radius) it declined until the torus's reversed meridian arcs
        // found their crossings with the box (PK-750: a reversed edge's
        // window errored and every crossing was dropped); it now serves
        // validator-clean with the same volume, every rim an exact arc as
        // at every finer tolerance (PK-779: circle recognition once
        // demanded a chord of 100 x the tolerance and left the 1e-2 rims
        // fitted).
        match popkern_bool::subtraction(torus, box_, t, None) {
            Ok(b) => {
                assert!(
                    b.validate_all().is_empty(),
                    "torus-box at {lin:e} must be validator-clean"
                );
                assert_eq!(
                    nurbs_edges(&b),
                    0,
                    "torus-box at {lin:e}: every rim is an exact arc"
                );
                let v = volume(&b);
                let rel = (v - truth_torus).abs() / truth_torus;
                assert!(
                    rel < 0.02,
                    "torus-box at {lin:e}: vol {v:.4} vs {truth_torus:.4}"
                );
            }
            Err(e) => panic!("torus-box at {lin:e} must serve, got decline: {e}"),
        }

        // The same torus cut IN its equator plane (z = 0): the tool's
        // bottom face contains the torus's own equator circles, so the
        // intersection records run along existing edges. The records
        // unify with those edges at the boolean tolerance and the kept
        // half torus welds to the annulus along them (PK-420: at 1e-3
        // an exact-node unification caught four of eight quarters and
        // the served body was two torus faces with curveless edges,
        // validator-clean, 4% of the volume). Serves through 1e-3
        // with the half-torus volume π²Rr².
        let torus = BRep::torus(1.0, 0.3, t).unwrap();
        let box_ = BRep::unit_cube(t)
            .unwrap()
            .transform(&Transform3::scale_uniform(3.0))
            .unwrap()
            .transform(&Transform3::translation(-1.5, -1.5, 0.0))
            .unwrap();
        match popkern_bool::subtraction(torus, box_, t, None) {
            Ok(b) => {
                // Serves clean at every tolerance in the envelope
                // (1e-2 included): a served body must be clean.
                assert!(
                    b.validate_all().is_empty(),
                    "torus-box(z=0) at {lin:e} must be validator-clean"
                );
                assert_eq!(
                    nurbs_edges(&b),
                    0,
                    "torus-box(z=0) at {lin:e}: every rim is an exact arc"
                );
                let v = volume(&b);
                let truth = std::f64::consts::PI.powi(2) * 0.09;
                let rel = (v - truth).abs() / truth;
                assert!(
                    rel < 0.02,
                    "torus-box(z=0) at {lin:e}: vol {v:.4} vs {truth:.4}"
                );
                let mesh = tessellate(&b, &TessellationOptions::default()).unwrap();
                assert_eq!(
                    open_edges(&mesh),
                    0,
                    "torus-box(z=0) at {lin:e}: watertight"
                );
            }
            Err(e) => {
                assert!(
                    lin >= 5e-3,
                    "torus-box(z=0) at {lin:e} must serve, got: {e}"
                );
            }
        }

        // Centered-sphere subtraction: serves clean through 1e-3; at
        // 1e-2 (tolerance = 10% of the cap height) the machinery
        // genuinely degrades and must DECLINE with the closure
        // diagnostics, never serve the torn body.
        let cube = BRep::unit_cube(t).unwrap();
        let sphere = BRep::sphere(0.6, t)
            .unwrap()
            .transform(&Transform3::translation(0.5, 0.5, 0.5))
            .unwrap();
        match popkern_bool::subtraction(cube, sphere, t, None) {
            Ok(b) => {
                assert!(
                    lin < 5e-3,
                    "sphere-cube at {lin:e} served — if the coarse regime now \
                     serves, it must be clean (below) and this pin updated"
                );
                assert!(
                    b.validate_all().is_empty(),
                    "sphere-cube at {lin:e} must be validator-clean"
                );
                let v = volume(&b);
                let rel = (v - truth_sphere).abs() / truth_sphere;
                assert!(
                    rel < 0.06,
                    "sphere-cube at {lin:e}: vol {v:.4} vs {truth_sphere:.4}"
                );
                // PK-418: the six cap rims (one split by the seam) and
                // the two seam pieces are exact arcs, and the mesh
                // above was taken with them in place — the sheet-
                // aligned hole rings keep the sphere face at its
                // 1e-4 triangle budget instead of over-filling.
                assert_eq!(
                    nurbs_edges(&b),
                    0,
                    "sphere-cube at {lin:e}: every rim and seam piece is an exact arc"
                );
            }
            Err(e) => {
                assert!(
                    lin >= 5e-3,
                    "sphere-cube at {lin:e} must serve, got decline: {e}"
                );
                let msg = format!("{e}");
                // The coarse regime's decline reason evolved: the
                // closure-diagnostic declines (loop gaps) were the
                // historical shape; the empty-serve lie gate
                // (PK-465) now catches the degradation EARLIER — the
                // garbage walk classifies everything away and the
                // boolean refuses the vanished body instead of
                // serving volume 0.0. Both are structured declines.
                assert!(
                    msg.contains("not closed")
                        || msg.contains("gap")
                        || msg.contains("deviates")
                        || msg.contains("empty result"),
                    "decline must carry structured diagnostics: {msg}"
                );
            }
        }
    }
}
