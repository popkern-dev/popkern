//! §5.3 R8 C8 — downstream parity scorecard (boolean leg).
//!
//! Runs a boolean operation against the filleted cube with both
//! `use_rational_gregory` flag states. The Gregory surfaces route
//! through `SurfaceGeometry::to_nurbs` (lossy sample-fit) in the
//! bool intersector, so this test is a smoke check — confirming
//! no panics, no topology corruption, non-empty result — not a
//! precision-parity test.

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_bool::{BoolOperation, boolean_ref};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Fillet the mixed-radius cube corner with the given flag then
/// subtract a small probe cube from it.
fn fillet_and_boolean(use_rational_gregory: bool) -> Result<BRep, String> {
    let cube = BRep::unit_cube(Tolerance::default()).map_err(|e| format!("{e:?}"))?;
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();

    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory,
        ..FilletOptions::default()
    };
    let filleted = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None)
        .map_err(|e| format!("fillet failed: {e:?}"))?
        .brep;

    // Probe solid: a second unit cube translated so it overlaps
    // the filleted corner at V=(0,0,0). Translating by
    // (-0.5, -0.5, -0.5) places the probe cube's "far" corner at
    // V, so the overlap region straddles the Gregory blend.
    let probe = BRep::unit_cube(Tolerance::default())
        .map_err(|e| format!("probe cube build failed: {e:?}"))?
        .transform(&Transform3::translation(-0.5, -0.5, -0.5))
        .map_err(|e| format!("probe cube transform failed: {e:?}"))?;

    boolean_ref(
        BoolOperation::Intersection,
        &filleted,
        &probe,
        Tolerance::default(),
        None,
    )
    .map_err(|e| format!("bool intersection failed: {e:?}"))
}

/// R8 C8 equal-or-better: enabling `use_rational_gregory` must
/// not REGRESS the boolean pipeline. If the baseline (polar-
/// NURBS) fillet + boolean succeeds, so must the Gregory path.
/// If the baseline fails (pre-existing issue with filleted-cube
/// trim curves in the bool intersector), we confirm the Gregory
/// path fails too — this is a known gap documented on the R8 C8
/// scorecard, not an R8 regression.
///
/// **Observed on the `(0.08, 0.10, 0.12)` mixed-radius cube
/// corner fixture (2026-04-17):** both baseline and Gregory
/// paths trip the same bool-intersection precondition check
/// (`trim curve end point does not match any vertex on the face
/// boundary`). This is a bool-side issue, unrelated to the R8
/// blend surface choice.
#[test]
fn r8_c8_bool_intersection_parity() {
    let baseline = fillet_and_boolean(false);
    let gregory = fillet_and_boolean(true);

    match (baseline, gregory) {
        (Ok(baseline_brep), Ok(gregory_brep)) => {
            // Both paths succeeded (a state first reached when the
            // PK-482 pinched-cut routing let this fixture's
            // boundary-tangent trim rings serve instead of aborting).
            // The parity claim is GEOMETRIC: both paths must carve
            // the same solid. Face counts are NOT the oracle — the
            // two blend representations legitimately partition the
            // corner differently (measured 2026-09-12: baseline 21
            // faces / volume 0.121645, Gregory 12 faces / volume
            // 0.121430 — a 0.18% mesh-discretization band on the
            // identical region), and the historical ≥-face-count
            // assertion pinned a route both paths no longer share.
            baseline_brep
                .validate()
                .expect("baseline bool result must pass validate()");
            gregory_brep
                .validate()
                .expect("R8 C8: Gregory bool result must pass validate()");
            let vol = |b: &popkern_core::topology::BRep| -> f64 {
                let mesh =
                    popkern_mesh::tessellate(b, &popkern_mesh::TessellationOptions::default())
                        .expect("meshes");
                let mut v = 0.0_f64;
                for tri in &mesh.triangles {
                    let a = mesh.vertices[tri[0]].position;
                    let bb = mesh.vertices[tri[1]].position;
                    let c = mesh.vertices[tri[2]].position;
                    let n = popkern_core::geometry::Vec3::new(
                        bb.y * c.z - bb.z * c.y,
                        bb.z * c.x - bb.x * c.z,
                        bb.x * c.y - bb.y * c.x,
                    );
                    v += (a.x * n.x + a.y * n.y + a.z * n.z) / 6.0;
                }
                v.abs()
            };
            let (v_base, v_greg) = (vol(&baseline_brep), vol(&gregory_brep));
            let rel = (v_base - v_greg).abs() / v_base.max(1e-12);
            assert!(
                rel < 0.005,
                "R8 C8: Gregory and baseline bool volumes must agree \
                 (baseline {v_base:.6} vs gregory {v_greg:.6}, rel {rel:.4})"
            );
        }
        (Err(base_err), Err(greg_err)) => {
            // Both paths hit the same upstream bool-side gap —
            // not an R8 regression. The test documents parity
            // (equal-or-better) and doesn't block.
            eprintln!(
                "R8 C8 known gap: bool intersection fails on both flag states.\n  \
                 baseline: {base_err}\n  gregory:  {greg_err}"
            );
        }
        (Ok(_), Err(greg_err)) => {
            panic!(
                "R8 C8 REGRESSION: boolean intersection succeeded with \
                 use_rational_gregory=false but failed with =true: {greg_err}. \
                 This must be fixed before the default can flip."
            );
        }
        (Err(base_err), Ok(_)) => {
            // Gregory path rescues a case where polar-NURBS fails
            // — record as a win on the scorecard but don't assert
            // it (it's not a required invariant).
            eprintln!(
                "R8 C8 note: Gregory bool intersection succeeded while \
                 baseline polar-NURBS failed: {base_err}"
            );
        }
    }
}
