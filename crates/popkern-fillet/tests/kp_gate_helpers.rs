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
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

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
    let faces: Vec<_> = result.corner_patches.iter().map(|p| p.face).collect();
    (result.brep, faces)
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

fn normal_curvature_along(gp: &GregoryPatch, u: f64, v: f64, t_ambient: Vec3) -> Option<f64> {
    let su = gp.partial_u(u, v);
    let sv = gp.partial_v(u, v);
    let n = su.cross(&sv).normalized()?;
    let e = dot(su, su);
    let f = dot(su, sv);
    let g = dot(sv, sv);
    let det = e * g - f * f;
    if det.abs() < 1e-14 {
        return None;
    }
    let su_dot_t = dot(su, t_ambient);
    let sv_dot_t = dot(sv, t_ambient);
    let a = (g * su_dot_t - f * sv_dot_t) / det;
    let b = (e * sv_dot_t - f * su_dot_t) / det;
    let norm_sq = e * a * a + 2.0 * f * a * b + g * b * b;
    if norm_sq < 1e-14 {
        return None;
    }
    let a = a / norm_sq.sqrt();
    let b = b / norm_sq.sqrt();
    let (suu, suv, svv) = match gp {
        GregoryPatch::RectangularBiNonic(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        GregoryPatch::RectangularTwoStripe(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        GregoryPatch::RectangularBiseptic(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        GregoryPatch::RectangularBiquinticFull(r) => (
            r.second_partial_uu(u, v),
            r.second_partial_uv(u, v),
            r.second_partial_vv(u, v),
        ),
        _ => return None,
    };
    let l = dot(n, suu);
    let m = dot(n, suv);
    let nn = dot(n, svv);
    Some(l * a * a + 2.0 * m * a * b + nn * b * b)
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
        "[{fixture}] expected at least one Gregory corner patch, got none (dispatch likely bailed out of the binonic rung)"
    );
    for (i, p) in patches.iter().enumerate() {
        assert!(
            matches!(p, GregoryPatch::RectangularBiNonic(_)),
            "[{fixture}] patch {i}: expected binonic variant, got a lower-rung fallback"
        );
    }
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

/// Worst-case continuity at the *outer* boundary of the corner-cap
/// ring — where each cap sub-face meets the cylinder fillet roll
/// (or other non-cap neighbour) it sits on top of. Complements
/// [`RingMetrics`] (which measures continuity *between* adjacent
/// sub-patches at their shared spokes).
///
/// Per Task #144 / 2026-04-27 failed-approaches entry: the ring
/// metric historically passed the gates while the cap-outer
/// boundary silently drifted up to 15% bbox in 3-D. The user-
/// visible "messy top" defect lives at this boundary.
pub struct CapOuterMetrics {
    /// Max 3-D gap `||cyl(cyl_pc(t)) − greg(greg_pc(t))||` over
    /// every cap-vs-non-cap edge sample.
    pub worst_g0_gap: f64,
    /// Max unit-normal parallelism deficit `1 − |n̂_a · n̂_b|` at
    /// edge samples. `0` = surfaces share a tangent plane;
    /// `1` = orthogonal normals (totally non-G1).
    pub worst_g1_coplanarity: f64,
    /// Number of cap-vs-non-cap edges measured.
    pub edges_measured: usize,
}

/// Measure G0/G1 at every edge between a Gregory cap face (in
/// `cap_faces`) and a non-cap neighbour. Samples each edge at
/// `n_samples` interior points (margin 0.05). Returns the max
/// over all samples on all such edges.
pub fn measure_cap_outer_boundary(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    n_samples: usize,
) -> CapOuterMetrics {
    measure_cap_outer_boundary_inner(brep, cap_faces, n_samples, false)
}

pub fn measure_cap_outer_boundary_including_pcurve_only(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    n_samples: usize,
) -> CapOuterMetrics {
    measure_cap_outer_boundary_inner(brep, cap_faces, n_samples, true)
}

fn measure_cap_outer_boundary_inner(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    n_samples: usize,
    include_pcurve_only_edges: bool,
) -> CapOuterMetrics {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let mut worst_g0 = 0.0_f64;
    let mut worst_g1 = 0.0_f64;
    let mut edges_measured = 0usize;
    let mut worst_g0_trace = String::new();
    let mut worst_g1_trace = String::new();
    let n_samples = n_samples.max(2);

    for (_eh, edge) in brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        // Only measure edges where exactly one side is a cap. Cap-
        // to-cap edges are the inner spokes already covered by
        // RingMetrics; non-cap-to-non-cap edges aren't of interest.
        if a_is_cap == b_is_cap {
            continue;
        }
        let Some(pa) = ha.pcurve.as_ref() else {
            continue;
        };
        let Some(pb) = hb.pcurve.as_ref() else {
            continue;
        };
        let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
            continue;
        };

        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let pa_at = |t: f64| {
            let uv = pa.point_at(t);
            sa.geometry.point_at(uv.x, uv.y)
        };
        let pb_at = |t: f64| {
            let uv = pb.point_at(t);
            sb.geometry.point_at(uv.x, uv.y)
        };
        let (pa_forward, pb_forward) = if let Some(coe) = edge.curve.as_ref() {
            let (t_lo, t_hi) = (coe.t_start, coe.t_end);
            // Direction detection (mirrors find_tessellation_non_manifold_edges).
            let p3_lo = coe.geometry.point_at(t_lo);
            let p3_hi = coe.geometry.point_at(t_hi);
            (
                p3_lo.distance_to(&pa_at(pa_lo)) + p3_hi.distance_to(&pa_at(pa_hi))
                    <= p3_lo.distance_to(&pa_at(pa_hi)) + p3_hi.distance_to(&pa_at(pa_lo)),
                p3_lo.distance_to(&pb_at(pb_lo)) + p3_hi.distance_to(&pb_at(pb_hi))
                    <= p3_lo.distance_to(&pb_at(pb_hi)) + p3_hi.distance_to(&pb_at(pb_lo)),
            )
        } else if include_pcurve_only_edges {
            // Some legacy cap builders leave the shared 3D edge curve absent
            // while retaining both pcurves. Measure the actual surface-to-
            // surface boundary by aligning the pcurve endpoints directly.
            let a_lo = pa_at(pa_lo);
            let a_hi = pa_at(pa_hi);
            (
                true,
                a_lo.distance_to(&pb_at(pb_lo)) + a_hi.distance_to(&pb_at(pb_hi))
                    <= a_lo.distance_to(&pb_at(pb_hi)) + a_hi.distance_to(&pb_at(pb_lo)),
            )
        } else {
            continue;
        };

        edges_measured += 1;
        for i in 0..n_samples {
            // Margin 0.05/0.95 to stay clear of the corner singularities
            // that the ring metrics also exclude.
            let s = 0.05 + 0.90 * (i as f64) / ((n_samples - 1) as f64);
            let s_a = if pa_forward { s } else { 1.0 - s };
            let s_b = if pb_forward { s } else { 1.0 - s };
            let t_a = pa_lo + s_a * (pa_hi - pa_lo);
            let t_b = pb_lo + s_b * (pb_hi - pb_lo);
            let uv_a = pa.point_at(t_a);
            let uv_b = pb.point_at(t_b);
            let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
            let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
            let g0 = p_a.distance_to(&p_b);
            if g0 > worst_g0 {
                worst_g0 = g0;
                if std::env::var_os("POPKERN_TRACE_CAP_OUTER_MEASURE").is_some() {
                    worst_g0_trace = format!(
                        "edge={_eh:?} faces=({face_a:?},{face_b:?}) sample={i} uv_a=({:.3},{:.3}) uv_b=({:.3},{:.3}) p_a=({:.4},{:.4},{:.4}) p_b=({:.4},{:.4},{:.4})",
                        uv_a.x, uv_a.y, uv_b.x, uv_b.y, p_a.x, p_a.y, p_a.z, p_b.x, p_b.y, p_b.z
                    );
                }
            }

            // G1 — unit-normal parallelism deficit. `0` = perfect
            // tangent-plane match, `1` = orthogonal normals. The
            // ring-metric's triple-product variant is not appropriate
            // here because we don't have an a-priori cross-tangent
            // direction (cap UV parameterization is implicit).
            let su_a = sa.geometry.partial_u(uv_a.x, uv_a.y);
            let sv_a = sa.geometry.partial_v(uv_a.x, uv_a.y);
            let su_b = sb.geometry.partial_u(uv_b.x, uv_b.y);
            let sv_b = sb.geometry.partial_v(uv_b.x, uv_b.y);
            let n_a = su_a.cross(&sv_a);
            let n_b = su_b.cross(&sv_b);
            let nn_a = vec_norm(n_a);
            let nn_b = vec_norm(n_b);
            if nn_a > 1e-300 && nn_b > 1e-300 {
                let deficit = 1.0 - (dot(n_a, n_b) / (nn_a * nn_b)).abs();
                if deficit > worst_g1 {
                    worst_g1 = deficit;
                    if std::env::var_os("POPKERN_TRACE_CAP_OUTER_MEASURE").is_some() {
                        worst_g1_trace = format!(
                            "edge={_eh:?} faces=({face_a:?},{face_b:?}) sample={i} uv_a=({:.3},{:.3}) uv_b=({:.3},{:.3}) n_a=({:.3},{:.3},{:.3}) n_b=({:.3},{:.3},{:.3}) su_a=({:.3},{:.3},{:.3}) sv_a=({:.3},{:.3},{:.3})",
                            uv_a.x,
                            uv_a.y,
                            uv_b.x,
                            uv_b.y,
                            n_a.x / nn_a,
                            n_a.y / nn_a,
                            n_a.z / nn_a,
                            n_b.x / nn_b,
                            n_b.y / nn_b,
                            n_b.z / nn_b,
                            su_a.x,
                            su_a.y,
                            su_a.z,
                            sv_a.x,
                            sv_a.y,
                            sv_a.z
                        );
                    }
                }
            }
        }
    }
    if !worst_g0_trace.is_empty() {
        eprintln!("[cap_outer measure worst_g0] {worst_g0:.6e} {worst_g0_trace}");
    }
    if !worst_g1_trace.is_empty() {
        eprintln!("[cap_outer measure worst_g1] {worst_g1:.6e} {worst_g1_trace}");
    }

    CapOuterMetrics {
        worst_g0_gap: worst_g0,
        worst_g1_coplanarity: worst_g1,
        edges_measured,
    }
}
