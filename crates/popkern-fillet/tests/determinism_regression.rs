//! §5.3 Phase-10a Issue B determinism regression gate.
//!
//! Runs the validation sweep's key scoring logic TWICE in the
//! same process and asserts the resulting metric vectors are
//! bit-identical. Any HashMap/HashSet iteration-order leak into
//! the output causes this test to fail, which is the definition
//! of non-deterministic behaviour we must eliminate.
//!
//! Uses a compact per-fixture metric: `(min_cos, mean_cos,
//! max_cos, max_position_gap)` for each variant rung. Compared
//! bit-for-bit (`f64::to_bits`) — tolerates nothing, because
//! genuine determinism means identical floating-point output.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct RawMetrics {
    min_cos_bits: u64,
    mean_cos_bits: u64,
    max_cos_bits: u64,
    max_gap_bits: u64,
    patch_count: usize,
}

fn fillet_cube_mixed_biseptic() -> RawMetrics {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    // Use BTreeSet for de-dup (edge appears via 2 half-edges) to
    // keep the incident-edges Vec stable across runs. A plain
    // HashSet here causes test-fixture non-determinism that
    // confounds pipeline-level determinism checks.
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    let radii = [0.08_f64, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None).unwrap();
    let patches: Vec<GregoryPatch> = result
        .corner_patches
        .iter()
        .filter_map(|cp| {
            result.brep.faces[cp.face].surface.as_ref().and_then(|sof| {
                if let SurfaceGeometry::Gregory(gp) = &sof.geometry {
                    Some(gp.clone())
                } else {
                    None
                }
            })
        })
        .collect();

    let mut min_cos: f64 = f64::INFINITY;
    let mut max_cos: f64 = f64::NEG_INFINITY;
    let mut sum_cos: f64 = 0.0;
    let mut count: usize = 0;
    let mut max_gap: f64 = 0.0;
    let n = patches.len();
    let ts: Vec<f64> = (1..=5).map(|k| (k as f64) / 6.0).collect();
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for &t in &ts {
            let t_b = 1.0 - t;
            let Some(na) = a.partial_u(1.0, t).cross(&a.partial_v(1.0, t)).normalized() else {
                continue;
            };
            let Some(nb) = b
                .partial_u(t_b, 0.0)
                .cross(&b.partial_v(t_b, 0.0))
                .normalized()
            else {
                continue;
            };
            let cos = (na.x * nb.x + na.y * nb.y + na.z * nb.z).abs();
            min_cos = min_cos.min(cos);
            max_cos = max_cos.max(cos);
            sum_cos += cos;
            count += 1;
            let pa = a.point_at(1.0, t);
            let pb = b.point_at(t_b, 0.0);
            let gap =
                ((pa.x - pb.x).powi(2) + (pa.y - pb.y).powi(2) + (pa.z - pb.z).powi(2)).sqrt();
            if gap > max_gap {
                max_gap = gap;
            }
        }
    }
    let mean_cos = if count > 0 {
        sum_cos / (count as f64)
    } else {
        0.0
    };
    RawMetrics {
        min_cos_bits: min_cos.to_bits(),
        mean_cos_bits: mean_cos.to_bits(),
        max_cos_bits: max_cos.to_bits(),
        max_gap_bits: max_gap.to_bits(),
        patch_count: patches.len(),
    }
}

/// Runs `fillet_cube_mixed_biseptic` twice, compares bit-for-bit.
///
/// This is the **definition** of the determinism bar we're
/// holding the pipeline to. If this passes, any two runs on the
/// same input produce identical output — the foundation for
/// reliable regression testing everywhere else.
///
/// Currently EXPECTED TO FAIL until the HashMap → BTreeMap
/// audit completes. When it passes, pin it to `#[test]` (not
/// `#[ignore]`) and the regression is locked down.
#[test]
fn bit_identical_output_on_repeat_run() {
    let a = fillet_cube_mixed_biseptic();
    let b = fillet_cube_mixed_biseptic();
    assert_eq!(
        a.patch_count, b.patch_count,
        "determinism failure: patch count differs ({} vs {})",
        a.patch_count, b.patch_count
    );
    assert_eq!(
        a.min_cos_bits,
        b.min_cos_bits,
        "determinism failure: min_cos differs (bits {:#x} vs {:#x}, values {:.17} vs {:.17})",
        a.min_cos_bits,
        b.min_cos_bits,
        f64::from_bits(a.min_cos_bits),
        f64::from_bits(b.min_cos_bits)
    );
    assert_eq!(
        a.mean_cos_bits,
        b.mean_cos_bits,
        "determinism failure: mean_cos differs ({:.17} vs {:.17})",
        f64::from_bits(a.mean_cos_bits),
        f64::from_bits(b.mean_cos_bits)
    );
    assert_eq!(
        a.max_cos_bits, b.max_cos_bits,
        "determinism failure: max_cos differs"
    );
    assert_eq!(
        a.max_gap_bits,
        b.max_gap_bits,
        "determinism failure: max_gap differs ({:.17} vs {:.17})",
        f64::from_bits(a.max_gap_bits),
        f64::from_bits(b.max_gap_bits)
    );
}
