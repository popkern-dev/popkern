//! Shared helpers for the Phase H.7 commercial-grade continuity
//! gate tests. Included into per-fixture test binaries via the
//! `#[path = "kp_gate_helpers.rs"] mod helpers;` idiom — keeps
//! the per-test files small enough to pass the 500-LoC structural
//! gate while centralising the metric definitions.
//!
//! Each test binary uses a subset of these helpers; the rest are
//! dead-code from the per-binary perspective. Top-level allow is
//! cleaner than per-item.
//!
//! See `kp_binonic_commercial_gates.rs` for the discussion of why
//! the triple-product coplanarity metric replaces the naive
//! unit-normal dot product.
#![allow(dead_code)]

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::GregoryPatch;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use popkern_fillet::{BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "support/kp_gate_cap_outer.rs"]
mod kp_gate_cap_outer;
#[allow(unused_imports)]
pub use kp_gate_cap_outer::{CapOuterMetrics, measure_cap_outer_boundary};
#[path = "support/kp_gate_shared_trims.rs"]
mod kp_gate_shared_trims;
#[allow(unused_imports)]
pub use kp_gate_shared_trims::measure_shared_gregory_trims_with_sources;
#[path = "support/kp_gate_curvature.rs"]
mod kp_gate_curvature;
#[allow(unused_imports)]
pub use kp_gate_curvature::normal_curvature_along;
#[path = "support/kp_gate_live_seam.rs"]
mod kp_gate_live_seam;
#[allow(unused_imports)]
pub use kp_gate_live_seam::{
    LiveSeamMetrics, measure_live_cap_seams, seam_cos_profile, seam_g2t_profile,
    welded_open_edge_count,
};

/// Commercial-grade gate thresholds. Discussion in
/// `kp_binonic_commercial_gates.rs` module header.
pub const N_SAMPLES: usize = 100;
pub const G0_GATE: f64 = 1e-10;
pub const G1_GATE: f64 = 5e-2;
pub const G2_GATE: f64 = 5e-2;

/// Deterministically-ordered incident edges at a vertex.
pub fn incident_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    brep.half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect()
}

/// First vertex with the requested valence.
pub fn pick_vertex_with_valence(brep: &BRep, target: usize) -> Option<VertexHandle> {
    brep.vertices
        .keys()
        .find(|&v| incident_edges(brep, v).len() == target)
}

/// Fillet a prepared `(brep, vertex, edges)` tuple with the full
/// deg5 → binonic promotion ladder. Returns `(brep, corner_faces)`.
pub fn fillet_corner(
    brep: &BRep,
    edges: &[EdgeHandle],
    radii: &[f64],
    use_kp_binonic: bool,
) -> (BRep, Vec<FaceHandle>) {
    fillet_corner_full(brep, edges, radii, use_kp_binonic, false)
}

/// Fillet a prepared tuple and collect cap faces by cap-surface semantics.
///
/// This is intentionally separate from [`fillet_corner`]: ring-continuity
/// tests need construction-order Gregory faces, while cap-outer tests need the
/// unordered set of actual cap surfaces, including PnS/sphere-cap faces that
/// are not reported as Gregory corner patches.
pub fn fillet_corner_cap_outer_faces(
    brep: &BRep,
    edges: &[EdgeHandle],
    radii: &[f64],
    use_kp_binonic: bool,
    use_kp_polar: bool,
) -> (BRep, Vec<FaceHandle>) {
    assert_eq!(edges.len(), radii.len());
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
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: use_kp_binonic,
        use_kp_polar_gregory: use_kp_polar,
        ..FilletOptions::default()
    };
    let result = fillet_edges(brep, &specs, &opts, Tolerance::default(), None).unwrap();
    let faces = cap_outer_faces_from_result(&result);
    (result.brep, faces)
}

/// Same as [`fillet_corner`], with an extra flag to enable the
/// Phase I polar K-P rung. `use_kp_polar = true` turns on the N ≥ 5
/// bi-6 → bi-9 promotion behind the `use_kp_polar_gregory` option.
pub fn fillet_corner_full(
    brep: &BRep,
    edges: &[EdgeHandle],
    radii: &[f64],
    use_kp_binonic: bool,
    use_kp_polar: bool,
) -> (BRep, Vec<FaceHandle>) {
    assert_eq!(edges.len(), radii.len());
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
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: use_kp_binonic,
        use_kp_polar_gregory: use_kp_polar,
        ..FilletOptions::default()
    };
    let result = fillet_edges(brep, &specs, &opts, Tolerance::default(), None).unwrap();
    let faces = cap_faces_from_result(&result);
    (result.brep, faces)
}

/// Cap faces in construction-report order (corner patches, then
/// exact-boundary committed faces), deduplicated without sorting.
/// Ring-continuity metrics are order-sensitive, so callers that
/// need adjacency must use this order, not SlotMap handle order.
pub fn cap_faces_from_result(result: &BlendResult) -> Vec<FaceHandle> {
    let mut faces = Vec::new();
    // Ring continuity metrics are order-sensitive. Preserve the
    // construction report order and deduplicate without sorting by
    // SlotMap handle, which is not a topological adjacency order.
    for face in result.corner_patches.iter().map(|patch| patch.face).chain(
        result
            .repair_report
            .apex_routes
            .iter()
            .flat_map(|route| route.exact_boundary_report.committed_faces.iter().copied()),
    ) {
        if !faces.contains(&face) {
            faces.push(face);
        }
    }
    faces
}

pub fn cap_outer_faces_from_result(result: &BlendResult) -> Vec<FaceHandle> {
    let mut faces = cap_faces_from_result(result);
    for (face, face_data) in result.brep.faces.iter() {
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if matches!(
            surface.geometry,
            SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Gregory(_)
        ) && !faces.contains(&face)
        {
            faces.push(face);
        }
    }
    faces
}

pub fn patch_of(brep: &BRep, face: FaceHandle) -> Option<GregoryPatch> {
    // Use `.get()` not index: at certain extreme-asymmetry fixtures
    // (cube_extreme_40 observed) the sliver-absorb or cap-insertion
    // passes in the fillet pipeline delete the corner face after
    // it's already been recorded in `result.corner_patches`, leaving
    // a stale handle. Panicking in the test helper obscures the
    // real (upstream) failure — silently skipping stale faces lets
    // the empty-corner-patches gate in `assert_all_binonic` fire
    // with its documented diagnostic.
    brep.faces.get(face)?.surface.as_ref().and_then(|sof| {
        if let SurfaceGeometry::Gregory(gp) = &sof.geometry {
            Some(gp.clone())
        } else {
            None
        }
    })
}

fn vec_norm(v: Vec3) -> f64 {
    (v.x * v.x + v.y * v.y + v.z * v.z).sqrt()
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

pub struct RingMetrics {
    pub worst_g0_gap: f64,
    pub worst_g1_coplanarity: f64,
    pub worst_g2_residual: f64,
}

/// Sample the ring on `t ∈ [0.05, 0.95]` at `n_samples` points.
/// Endpoint margins exclude the Gregory corner-singularity
/// transition zones — not a Phase H construction concern.
pub fn measure_ring(patches: &[GregoryPatch], n_samples: usize) -> RingMetrics {
    let n = patches.len();
    let mut worst_g0 = 0.0_f64;
    let mut worst_g1 = 0.0_f64;
    let mut worst_g2 = 0.0_f64;
    let ts: Vec<f64> = (0..n_samples)
        .map(|k| 0.05 + 0.90 * (k as f64) / ((n_samples - 1) as f64))
        .collect();
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for &t in &ts {
            let t_b = 1.0 - t;

            let pa = a.point_at(1.0, t);
            let pb = b.point_at(t_b, 0.0);
            let gap =
                ((pa.x - pb.x).powi(2) + (pa.y - pb.y).powi(2) + (pa.z - pb.z).powi(2)).sqrt();
            if gap > worst_g0 {
                worst_g0 = gap;
            }

            let xtan_a = a.partial_u(1.0, t);
            let xtan_b = b.partial_v(t_b, 0.0);
            let t_spoke = a.partial_v(1.0, t);
            let numer = dot(xtan_a.cross(&xtan_b), t_spoke).abs();
            let denom = vec_norm(xtan_a) * vec_norm(xtan_b) * vec_norm(t_spoke);
            if denom > 1e-300 {
                let coplan = numer / denom;
                if coplan > worst_g1 {
                    worst_g1 = coplan;
                }
            }

            let Some(t_unit) = t_spoke.normalized() else {
                continue;
            };
            let Some(ka) = normal_curvature_along(a, 1.0, t, t_unit) else {
                continue;
            };
            let Some(kb) = normal_curvature_along(b, t_b, 0.0, t_unit) else {
                continue;
            };
            // Fold the two charts' natural normals to a consistent
            // side before comparing II: a chart-handedness reversal
            // flips the SIGN of the measured normal curvature while
            // the geometry stays curvature-continuous (measured on
            // the PK-240 balanced rings: kb = -ka to 2-3 digits with
            // |dk| under 4e-2, while the unfolded comparison read
            // 1.2e0 — the same unfolded-normal artifact class as
            // PK-244's flags). Chart-reversal bands themselves are
            // recorded separately; this metric measures GEOMETRY.
            let na = a.partial_u(1.0, t).cross(&a.partial_v(1.0, t));
            let nb = b.partial_u(t_b, 0.0).cross(&b.partial_v(t_b, 0.0));
            // Conditioning guard: II is unreadable through a chart
            // whose basis is near-collinear (measured on the PK-240
            // balanced rings: the worst "G2" sample sat at a 0.0-deg
            // basis angle with healthy factor magnitudes — the
            // absolute det floor cannot catch it). Skip, don't clamp.
            let sin_a = vec_norm(na)
                / (vec_norm(a.partial_u(1.0, t)) * vec_norm(a.partial_v(1.0, t))).max(1e-300);
            let sin_b = vec_norm(nb)
                / (vec_norm(b.partial_u(t_b, 0.0)) * vec_norm(b.partial_v(t_b, 0.0))).max(1e-300);
            if sin_a < 5.0e-2 || sin_b < 5.0e-2 {
                continue;
            }
            let denom = vec_norm(na) * vec_norm(nb);
            if denom < 1e-300 {
                continue;
            }
            let cos = dot(na, nb) / denom;
            // Convergence guard: near a chart-degeneracy crossing the
            // natural normal spins through zero and II reads noise —
            // the G1 gate (coplanarity, orientation-free) separately
            // covers plane agreement, so a low |cos| here can only be
            // chart noise, never geometry. Skip, don't clamp.
            if cos.abs() < 0.9 {
                continue;
            }
            let kb = if cos < 0.0 { -kb } else { kb };
            let kmax = ka.abs().max(kb.abs());
            let resid = (ka - kb).abs() / (1.0 + kmax);
            if resid > worst_g2 {
                worst_g2 = resid;
            }
        }
    }
    RingMetrics {
        worst_g0_gap: worst_g0,
        worst_g1_coplanarity: worst_g1,
        worst_g2_residual: worst_g2,
    }
}

/// Collect Gregory sub-patches from the corner faces of a filleted
/// BRep. Panics if any face lacks a Gregory surface — every fixture
/// in the matrix is known to produce Gregory corners.
pub fn collect_patches(brep: &BRep, faces: &[FaceHandle]) -> Vec<GregoryPatch> {
    faces.iter().filter_map(|&f| patch_of(brep, f)).collect()
}

/// Assert every patch is the `RectangularBiNonic` variant —
/// proving the K-P promotion fired on this fixture. Also requires
/// at least one patch so an empty corner-patch list (dispatch
/// didn't produce any Gregory geometry — a common failure mode
/// on unsupported valences) is caught instead of passing the
/// `for` loop vacuously.
pub fn assert_all_binonic(patches: &[GregoryPatch], fixture: &str) {
    assert!(
        !patches.is_empty(),
        "[{fixture}] expected at least one Gregory corner patch, got none (dispatch likely bailed out before a measurable K-P ring)"
    );
    for (i, p) in patches.iter().enumerate() {
        assert!(
            matches!(p, GregoryPatch::RectangularBiNonic(_)),
            "[{fixture}] patch {i}: expected binonic variant, got a lower-rung fallback"
        );
    }
}

pub fn all_binonic(patches: &[GregoryPatch]) -> bool {
    !patches.is_empty()
        && patches
            .iter()
            .all(|patch| matches!(patch, GregoryPatch::RectangularBiNonic(_)))
}

/// Assert all three commercial gates (G0, G1, G2) pass on the
/// given metrics. Prints the metrics for logs.
pub fn assert_commercial_gates(fixture: &str, m: &RingMetrics) {
    println!(
        "binonic {fixture}: g0={:.3e} g1={:.3e} g2={:.3e}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.worst_g2_residual
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "[{fixture}] binonic G0 gap {:.3e} exceeds {:.1e} gate",
        m.worst_g0_gap,
        G0_GATE
    );
    assert!(
        m.worst_g1_coplanarity < G1_GATE,
        "[{fixture}] binonic G1 coplanarity {:.3e} exceeds {:.1e} gate (commercial)",
        m.worst_g1_coplanarity,
        G1_GATE
    );
    assert!(
        m.worst_g2_residual < G2_GATE,
        "[{fixture}] binonic G2 residual {:.3e} exceeds {:.1e} gate (commercial)",
        m.worst_g2_residual,
        G2_GATE
    );
}

/// Assert the commercial gates (G0, G1, G2) over the LIVE seam
/// topology — the layout-honest sibling of [`assert_commercial_gates`]
/// for multi-patch cap commits (the setback fan), where the legacy
/// ring convention misreads the frame (PK-087: 5.2e-2 on a fan whose
/// true seams sit at machine zero). Same thresholds: G0 is side-vs-
/// shared-curve (the kernel trim contract), G1 the normal cross-
/// product norm, G2 the normal-curvature jump along the seam.
pub fn assert_commercial_gates_live(fixture: &str, m: &LiveSeamMetrics) {
    println!(
        "live-seam {fixture}: seams={} g0={:.3e} g1={:.3e} g2={:.3e} g2t={:.3e}@f={:.2} ({:?})",
        m.seam_count,
        m.worst_g0_gap,
        m.worst_g1_normal,
        m.worst_g2_residual,
        m.worst_g2_transverse,
        m.worst_g2t_fraction,
        m.worst_g2t_edge
    );
    assert!(
        m.seam_count > 0,
        "[{fixture}] expected live cap-cap seams to measure, found none"
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "[{fixture}] live-seam G0 gap {:.3e} exceeds {:.1e} gate",
        m.worst_g0_gap,
        G0_GATE
    );
    assert!(
        m.worst_g1_normal < G1_GATE,
        "[{fixture}] live-seam G1 normal misalignment {:.3e} exceeds {:.1e} gate (commercial)",
        m.worst_g1_normal,
        G1_GATE
    );
    assert!(
        m.worst_g2_residual < G2_GATE,
        "[{fixture}] live-seam G2 residual {:.3e} exceeds {:.1e} gate (commercial)",
        m.worst_g2_residual,
        G2_GATE
    );
    // The transverse term is the honest curvature-continuity gate —
    // the along-seam residual above is near-vacuous once G0+G1 hold
    // (Meusnier; see LiveSeamMetrics docs). Before the PK-187
    // rim-solve G2-row pinning this measured 4.3e-1..1.7e0 on every
    // polar fixture while the along-seam gate read 1e-5 "commercial".
    assert!(
        m.worst_g2_transverse < G2_GATE,
        "[{fixture}] live-seam TRANSVERSE G2 residual {:.3e} exceeds {:.1e} gate \
         (commercial) at f={:.2} on {:?}",
        m.worst_g2_transverse,
        G2_GATE,
        m.worst_g2t_fraction,
        m.worst_g2t_edge
    );
}

/// Assert binonic is within 2× of two-stripe on every metric
/// (floor at 1e-8 / 1e-12 to absorb floating-point noise when
/// both values are essentially machine-ε).
pub fn assert_no_regression(fixture: &str, mt: &RingMetrics, mb: &RingMetrics) {
    println!(
        "two_stripe {fixture}: g0={:.3e} g1={:.3e} g2={:.3e}",
        mt.worst_g0_gap, mt.worst_g1_coplanarity, mt.worst_g2_residual
    );
    println!(
        "binonic    {fixture}: g0={:.3e} g1={:.3e} g2={:.3e}",
        mb.worst_g0_gap, mb.worst_g1_coplanarity, mb.worst_g2_residual
    );
    assert!(
        mb.worst_g0_gap <= 2.0 * mt.worst_g0_gap.max(1e-12),
        "[{fixture}] binonic G0 regression: {:.3e} vs {:.3e}",
        mb.worst_g0_gap,
        mt.worst_g0_gap
    );
    assert!(
        mb.worst_g1_coplanarity <= 2.0 * mt.worst_g1_coplanarity.max(1e-8),
        "[{fixture}] binonic G1 regression: {:.3e} vs {:.3e}",
        mb.worst_g1_coplanarity,
        mt.worst_g1_coplanarity
    );
    assert!(
        mb.worst_g2_residual <= 2.0 * mt.worst_g2_residual.max(1e-8),
        "[{fixture}] binonic G2 regression: {:.3e} vs {:.3e}",
        mb.worst_g2_residual,
        mt.worst_g2_residual
    );
}
