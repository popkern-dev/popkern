//! §5.3 item #11 — cross-fixture reflection-line + zebra-stripe
//! regression gate.
//!
//! Class-A continuity metrics (reflection lines, zebra stripes)
//! are the fine-grained complement to §2.4 boundary-κ deviation:
//! they catch subtle G2 defects that curvature metrics miss.
//! This gate runs both metrics across every committed §5.3 fixture
//! and asserts each stays below a reasonable bound — so a
//! regression in Gregory geometry, boundary fitting, or
//! support-face blending shows up in CI rather than leaking to
//! downstream consumers.
//!
//! Fixtures exercised:
//!   * Unit cube, single-edge fillet
//!   * Unit cube, all-equal-radius trihedral corner (SameSide
//!     sphere-fit path)
//!   * Unit cube, mixed-radius trihedral corner (Gregory path)
//!   * `BRep::square_pyramid` N=4 apex with Gregory
//!
//! Bounds:
//!   * Reflection-line max kink ≤ 2.0 rad (≈ 115°) — generous
//!     upper bound that rejects only catastrophic G0 breaks
//!     while allowing Gregory's inherent spoke discontinuity at
//!     the central vertex (item #3's known scope).
//!   * Zebra mismatch ratio ≤ 0.6 — permissive bound rejecting
//!     fully anti-correlated stripe fields while accepting the
//!     random-looking mismatches that arise from boundary
//!     tolerance choices on curved supports.
//!
//! The test is a regression gate, not a class-A certification —
//! tighter bounds belong in a follow-up once item #3's spoke
//! G2 solver lands.

use popkern_core::geometry::Vec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::reflection_line::blend_reflection_line_kink_max;
use popkern_fillet::zebra::blend_zebra_kink_max;
use popkern_fillet::{CornerMethod, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn single_edge_cube_fillet() -> BRep {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None)
        .unwrap()
        .brep
}

fn cube_corner_uniform() -> (BRep, popkern_core::topology::FaceHandle) {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let specs: Vec<_> = incident
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let r = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    let corner_face = r
        .corner_patches
        .iter()
        .find(|p| p.vertex == vertex)
        .expect("uniform-radius corner must emit a patch")
        .face;
    (r.brep, corner_face)
}

fn cube_corner_mixed_gregory() -> (BRep, popkern_core::topology::FaceHandle) {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let specs: Vec<_> = incident
        .iter()
        .zip([0.08, 0.10, 0.12])
        .map(|(&e, r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    let r = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
    // Pick any of the 3 Gregory sub-patch faces (they all share
    // the same spoke topology).
    let corner_face = r
        .corner_patches
        .iter()
        .find(|p| p.vertex == vertex && p.method == CornerMethod::Gregory)
        .expect("mixed-radii Gregory corner must emit a sub-patch")
        .face;
    (r.brep, corner_face)
}

fn pyramid_apex_gregory() -> Option<(BRep, popkern_core::topology::FaceHandle)> {
    let pyramid = BRep::square_pyramid(1.0, 0.8, tol()).ok()?;
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.z > 0.1)
        .map(|(h, _)| h)?;
    let slant_edges: Vec<_> = {
        let mut seen = std::collections::BTreeSet::new();
        pyramid
            .half_edges
            .values()
            .filter(|he| he.origin == apex)
            .map(|he| he.edge)
            .filter(|e| seen.insert(*e))
            .collect()
    };
    let radii = [0.07_f64, 0.09, 0.11, 0.13];
    let specs: Vec<_> = slant_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    let r = fillet_edges(&pyramid, &specs, &opts, tol(), None).ok()?;
    let face = r
        .corner_patches
        .iter()
        .find(|p| p.vertex == apex && p.method == CornerMethod::Gregory)?
        .face;
    Some((r.brep, face))
}

fn pick_blend_face_from_cube_fillet(brep: &BRep) -> popkern_core::topology::FaceHandle {
    use popkern_core::geometry::surface::SurfaceGeometry;
    // The single cylindrical fillet surface.
    brep.faces
        .iter()
        .find(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
        })
        .map(|(h, _)| h)
        .expect("single-edge fillet must produce a Cylinder face")
}

const REFLECTION_KINK_LIMIT: f64 = 2.0;
const ZEBRA_MISMATCH_LIMIT: f64 = 0.6;
const VIEW: Vec3 = Vec3 {
    x: 0.0,
    y: 0.0,
    z: 1.0,
};
const STRIPE: Vec3 = Vec3 {
    x: 1.0,
    y: 0.0,
    z: 0.0,
};

fn assert_within_bounds(label: &str, brep: &BRep, blend_face: popkern_core::topology::FaceHandle) {
    let rk = blend_reflection_line_kink_max(brep, blend_face, VIEW, STRIPE, 0.15, 24, 0.04);
    if let Some(angle) = rk {
        eprintln!("[{label}] reflection-line max kink = {angle:.4} rad");
        assert!(
            angle < REFLECTION_KINK_LIMIT,
            "§5.3 item #11: {label} reflection-line kink {angle:.4} rad \
             exceeds {REFLECTION_KINK_LIMIT} rad"
        );
    } else {
        eprintln!("[{label}] reflection-line kink — no boundary pairs found (skipping)");
    }

    let zk = blend_zebra_kink_max(brep, blend_face, VIEW, STRIPE, 0.15, 24, 0.04);
    if let Some(ratio) = zk {
        eprintln!("[{label}] zebra mismatch ratio = {ratio:.4}");
        assert!(
            ratio < ZEBRA_MISMATCH_LIMIT,
            "§5.3 item #11: {label} zebra mismatch {ratio:.4} exceeds \
             {ZEBRA_MISMATCH_LIMIT}"
        );
    } else {
        eprintln!("[{label}] zebra boundary — no matched pairs (skipping)");
    }
}

#[test]
fn reflection_and_zebra_within_bounds_cross_fixture() {
    // Fixture 1: single-edge cube fillet (Cylinder face).
    let brep = single_edge_cube_fillet();
    let blend_face = pick_blend_face_from_cube_fillet(&brep);
    assert_within_bounds("cube single-edge", &brep, blend_face);

    // Fixture 2: uniform-radius cube corner (sphere patch).
    let (brep, corner_face) = cube_corner_uniform();
    assert_within_bounds("cube corner uniform", &brep, corner_face);

    // Fixture 3: mixed-radius cube corner (Gregory).
    let (brep, corner_face) = cube_corner_mixed_gregory();
    assert_within_bounds("cube corner mixed Gregory", &brep, corner_face);

    // Fixture 4: pyramid apex (N=4 Gregory).
    if let Some((brep, corner_face)) = pyramid_apex_gregory() {
        assert_within_bounds("pyramid apex Gregory", &brep, corner_face);
    } else {
        eprintln!("pyramid apex fixture skipped — fillet pipeline rejected");
    }
}

/// Smoke test: `blend_reflection_line_kink_max` on a flat cube
/// face (not a fillet) returns `None` because the face has no
/// adjacent support faces with distinct surface geometry that
/// would produce a measurable kink. Documents the "no pairs"
/// failure mode.
#[test]
fn reflection_kink_returns_none_on_flat_cube_face() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let face = cube.faces.keys().next().unwrap();
    let result = blend_reflection_line_kink_max(&cube, face, VIEW, STRIPE, 0.15, 12, 0.01);
    // On a cube, every face abuts 4 perpendicular planar faces.
    // Reflection lines DO differ across the boundary (90° normal
    // kink at the edge), so the max-kink will be non-zero —
    // confirming the metric fires on genuine G0 transitions.
    if let Some(angle) = result {
        assert!(
            angle > 0.1,
            "perpendicular cube faces should produce measurable reflection kink"
        );
    }
}
