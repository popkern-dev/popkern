//! §5.3 Phase-10a Issue B determinism regression gate (PK-178).
//!
//! Runs the fillet pipeline TWICE in the same process and asserts
//! the outputs are bit-identical. Any HashMap/HashSet
//! iteration-order leak into the output causes this test to fail,
//! which is the definition of non-deterministic behaviour we must
//! eliminate.
//!
//! Two layers of comparison, both tolerating nothing:
//! - a compact Gregory corner-patch metric `(min_cos, mean_cos,
//!   max_cos, max_position_gap)` compared via `f64::to_bits`, and
//! - a whole-body fingerprint hashing the serialized output
//!   `BRep` — topology arenas, handles, and every geometry f64 —
//!   so routes that commit no `corner_patches` records (e.g. the
//!   direct-apex partition, PK-099) are still fully pinned.
//!
//! Fixtures cover two apex families through two dispatch routes:
//! cube N=3 (biseptic Gregory rungs) and square-pyramid N=4 mixed
//! (default commercial options → direct-apex partition), so
//! determinism is not asserted on a single incidentally-stable
//! fixture.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

/// FNV-1a over the serde_json serialization of the output body.
///
/// serde_json prints f64 via shortest-roundtrip (ryu), which is
/// injective on finite values, so two bodies hash equal only if
/// every coordinate, control point, and handle matches exactly.
fn body_fingerprint(brep: &BRep) -> u64 {
    let bytes = serde_json::to_vec(brep).expect("output BRep must serialize");
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        hash ^= u64::from(b);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// Deterministically ordered incident edges at `vertex` (BTreeSet
/// dedup — a HashSet here would inject fixture-level variance that
/// confounds pipeline-level determinism checks).
fn incident_edges(brep: &BRep, vertex: VertexHandle) -> Vec<popkern_core::topology::EdgeHandle> {
    brep.half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct RawMetrics {
    min_cos_bits: u64,
    mean_cos_bits: u64,
    max_cos_bits: u64,
    max_gap_bits: u64,
    patch_count: usize,
    body_hash: u64,
}

fn fillet_cube_mixed_biseptic() -> RawMetrics {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let radii = [0.08_f64, 0.10, 0.12];
    let specs: Vec<FilletSpec> = edges
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
        body_hash: body_fingerprint(&result.brep),
    }
}

/// Square-pyramid N=4 apex, mixed radii, default commercial
/// options — dispatches through the direct-apex partition
/// (PK-099), a different route from the cube's biseptic Gregory
/// rungs. Partition caps commit no `corner_patches` records, so
/// this fixture is pinned by the whole-body fingerprint.
fn fillet_pyramid4_mixed_default() -> u64 {
    let pyramid = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let apex = pyramid
        .vertices
        .keys()
        .find(|&v| incident_edges(&pyramid, v).len() == 4)
        .expect("square_pyramid has a valence-4 apex");
    let edges = incident_edges(&pyramid, apex);
    let radii = [0.06_f64, 0.08, 0.10, 0.08];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions::default();
    let result = fillet_edges(&pyramid, &specs, &opts, Tolerance::default(), None).unwrap();
    assert!(
        result.brep.faces.len() > pyramid.faces.len(),
        "pyramid fillet must actually commit blend faces (got {} vs input {})",
        result.brep.faces.len(),
        pyramid.faces.len()
    );
    body_fingerprint(&result.brep)
}

/// Runs `fillet_cube_mixed_biseptic` twice, compares bit-for-bit.
///
/// This is the **definition** of the determinism bar we're
/// holding the pipeline to. If this passes, any two runs on the
/// same input produce identical output — the foundation for
/// reliable regression testing everywhere else.
///
/// Pinned as a default `#[test]` since PK-178: the HashMap →
/// BTreeMap audit landed and the gate passes; it now locks the
/// regression down.
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
    assert_eq!(
        a.body_hash, b.body_hash,
        "determinism failure: serialized output body differs (hash {:#x} vs {:#x})",
        a.body_hash, b.body_hash
    );
}

/// Same bar on the square-pyramid N=4 mixed fixture through the
/// default commercial route (direct-apex partition), so pipeline
/// determinism is not asserted on a single incidentally-stable
/// fixture or route (PK-178 acceptance).
#[test]
fn bit_identical_output_on_repeat_run_pyramid4_mixed() {
    let a = fillet_pyramid4_mixed_default();
    let b = fillet_pyramid4_mixed_default();
    assert_eq!(
        a, b,
        "determinism failure: pyramid4_mixed serialized output body differs (hash {a:#x} vs {b:#x})"
    );
}

/// The 12-edge G2 cube fillet — the PK-278 two-phase corner path:
/// every corner's KP20 cap solve is precomputed on the parallel
/// presolve and consumed by the serial apply. Whole-body
/// fingerprint of the served result.
fn fillet_cube_g2_ladder() -> u64 {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .take(12)
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let opts = FilletOptions {
        continuity: popkern_fillet::ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None)
        .expect("12-edge G2 cube fillet serves");
    assert!(
        result.brep.faces.len() > cube.faces.len(),
        "G2 cube fillet must commit blend faces"
    );
    body_fingerprint(&result.brep)
}

/// PK-278 acceptance (row 416): the parallel corner presolve and
/// the parallel tail measurement passes must not perturb a single
/// output bit run-to-run.
#[test]
fn bit_identical_output_on_repeat_run_g2_cube() {
    let a = fillet_cube_g2_ladder();
    let b = fillet_cube_g2_ladder();
    assert_eq!(
        a, b,
        "determinism failure: G2 cube serialized output body differs (hash {a:#x} vs {b:#x})"
    );
}

/// PK-278 acceptance (row 416): the output must be bit-identical
/// at any worker-thread count. The presolve is a pure solve gated
/// on exact input equality and every parallel measurement pass
/// uses an indexed ordered collect, so 1 thread and 8 threads must
/// reproduce the serial sequence exactly.
#[test]
fn thread_count_does_not_change_g2_cube() {
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("1-thread pool builds")
        .install(fillet_cube_g2_ladder);
    let eight = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .expect("8-thread pool builds")
        .install(fillet_cube_g2_ladder);
    assert_eq!(
        one, eight,
        "determinism failure: G2 cube output differs between 1-thread and 8-thread pools \
         (hash {one:#x} vs {eight:#x})"
    );
}
