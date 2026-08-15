//! Gate metrics for the chart-orientation repair (PK-249) —
//! split from `three_quad_chart_orientation` for the 500-line
//! structural cap. Owns the orientation-independent measurement
//! helpers: fold detection, characteristic length, cap-vs-fillet
//! outer tangency, and the honest (handedness-folded,
//! singular-skipping) ring G2.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;

pub(super) fn any_self_intersects(ring: &[RectangularGregoryBiNonic], tolerance: f64) -> bool {
    ring.iter().any(|patch| {
        GregoryPatch::RectangularBiNonic(Box::new(patch.clone()))
            .self_intersects_with_tolerance(tolerance)
    })
}

pub(super) fn characteristic_length(ring: &[RectangularGregoryBiNonic]) -> f64 {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for patch in ring {
        for row in &patch.grid {
            for p in row {
                for (k, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                    min[k] = min[k].min(v);
                    max[k] = max[k].max(v);
                }
            }
        }
    }
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    (dx * dx + dy * dy + dz * dz).sqrt().max(1e-9)
}

/// Cap-vs-fillet G1 at one `u=0` boundary sample: the sine of the
/// angle between the cap normal and the projected fillet normal.
pub(super) fn outer_fillet_sin(
    binonic: &RectangularGregoryBiNonic,
    fillet: &SurfaceGeometry,
    t: f64,
) -> Option<f64> {
    let n_cap = binonic.partial_u(0.0, t).cross(&binonic.partial_v(0.0, t));
    let p = binonic.point_at(0.0, t);
    let (u_f, v_f) = project_point_to_uv(fillet, &p)?;
    let n_f = fillet
        .partial_u(u_f, v_f)
        .cross(&fillet.partial_v(u_f, v_f));
    let denom = n_cap.length() * n_f.length();
    (denom > 1e-300).then(|| n_cap.cross(&n_f).length() / denom)
}

pub(super) fn worst_outer_fillet_sin(
    ring: &[RectangularGregoryBiNonic],
    fillets: &[SurfaceGeometry],
) -> f64 {
    let mut worst = 0.0_f64;
    for (k, binonic) in ring.iter().enumerate() {
        let Some(fillet) = fillets.get(k) else {
            continue;
        };
        for s in 0..25 {
            let t = 0.05 + 0.90 * (s as f64) / 24.0;
            if let Some(sin) = outer_fillet_sin(binonic, fillet, t) {
                worst = worst.max(sin);
            }
        }
    }
    worst
}

/// Honest ring G2: normal curvature agreement across shared spokes
/// with chart handedness folded to a consistent side and
/// singular-basis / non-converged samples skipped — the geometry
/// metric, matching the commercial-gate tests' convention.
pub(super) fn worst_ring_folded_g2(ring: &[RectangularGregoryBiNonic]) -> f64 {
    let n = ring.len();
    let mut worst = 0.0_f64;
    for k in 0..n {
        let a = &ring[k];
        let b = &ring[(k + 1) % n];
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            if let Some(resid) = g2_residual_at(a, b, t) {
                worst = worst.max(resid);
            }
        }
    }
    worst
}

/// The derivative stack one folded-G2 sample reads: first and
/// second partials of both patches at their matched spoke points.
/// The bi-nonic patch is a plain tensor Bézier, so every entry is
/// LINEAR in the control points — the polish Jacobian exploits this
/// by perturbing the STACK by exact basis sensitivities instead of
/// re-evaluating the surfaces per parameter column (PK-282).
#[derive(Clone, Copy)]
pub(super) struct G2SampleStack {
    /// Patch `a` at `(1, t)`: `[su, sv, suu, suv, svv]`.
    pub a: [Vec3; 5],
    /// Patch `b` at `(1−t, 0)`: `[su, sv, suu, suv, svv]`.
    pub b: [Vec3; 5],
}

/// Gather the folded-G2 derivative stack for the spoke shared by
/// patch `a` (its `u=1` side at `t`) and patch `b` (its `v=0` side
/// at `1−t`).
pub(super) fn g2_sample_stack(
    a: &RectangularGregoryBiNonic,
    b: &RectangularGregoryBiNonic,
    t: f64,
) -> G2SampleStack {
    let t_b = 1.0 - t;
    G2SampleStack {
        a: [
            a.partial_u(1.0, t),
            a.partial_v(1.0, t),
            a.second_partial_uu(1.0, t),
            a.second_partial_uv(1.0, t),
            a.second_partial_vv(1.0, t),
        ],
        b: [
            b.partial_u(t_b, 0.0),
            b.partial_v(t_b, 0.0),
            b.second_partial_uu(t_b, 0.0),
            b.second_partial_uv(t_b, 0.0),
            b.second_partial_vv(t_b, 0.0),
        ],
    }
}

/// One folded normal-curvature agreement sample across the spoke
/// shared by patch `a` (its `u=1` side at `t`) and patch `b` (its
/// `v=0` side at `1−t`). `None` when the sample is unreadable
/// (singular basis or non-converged normals) — skipped, matching
/// the honest-G2 convention.
pub(super) fn g2_residual_at(
    a: &RectangularGregoryBiNonic,
    b: &RectangularGregoryBiNonic,
    t: f64,
) -> Option<f64> {
    g2_residual_from_stack(&g2_sample_stack(a, b, t))
}

/// The pure folded-G2 formula over a gathered derivative stack —
/// branch gates included, no surface evaluations.
pub(super) fn g2_residual_from_stack(stack: &G2SampleStack) -> Option<f64> {
    let [su_a, sv_a, ..] = stack.a;
    let [su_b, sv_b, ..] = stack.b;
    let t_unit = sv_a.normalized()?;
    let ka = normal_curvature_from_stack(&stack.a, t_unit)?;
    let kb = normal_curvature_from_stack(&stack.b, t_unit)?;
    let na = su_a.cross(&sv_a);
    let nb = su_b.cross(&sv_b);
    let sin_a = na.length() / (su_a.length() * sv_a.length()).max(1e-300);
    let sin_b = nb.length() / (su_b.length() * sv_b.length()).max(1e-300);
    if sin_a < 5.0e-2 || sin_b < 5.0e-2 {
        return None;
    }
    let denom = na.length() * nb.length();
    if denom < 1e-300 {
        return None;
    }
    let cos = na.dot(&nb) / denom;
    if cos.abs() < 0.9 {
        return None;
    }
    let kb = if cos < 0.0 { -kb } else { kb };
    let kmax = ka.abs().max(kb.abs());
    Some((ka - kb).abs() / (1.0 + kmax))
}

/// Normal curvature along the ambient direction `t_ambient`
/// (projected into the tangent plane and normalized through the
/// first fundamental form), over a `[su, sv, suu, suv, svv]` stack.
fn normal_curvature_from_stack(d: &[Vec3; 5], t_ambient: Vec3) -> Option<f64> {
    let [su, sv, suu, suv, svv] = *d;
    let normal = su.cross(&sv).normalized()?;
    let e = su.dot(&su);
    let f = su.dot(&sv);
    let g = sv.dot(&sv);
    let det = e * g - f * f;
    if det.abs() < 1e-14 {
        return None;
    }
    let su_dot_t = su.dot(&t_ambient);
    let sv_dot_t = sv.dot(&t_ambient);
    let a = (g * su_dot_t - f * sv_dot_t) / det;
    let b = (e * sv_dot_t - f * su_dot_t) / det;
    let norm_sq = e * a * a + 2.0 * f * a * b + g * b * b;
    if norm_sq < 1e-14 {
        return None;
    }
    let a = a / norm_sq.sqrt();
    let b = b / norm_sq.sqrt();
    let l = suu.dot(&normal);
    let m = suv.dot(&normal);
    let nn = svv.dot(&normal);
    Some(l * a * a + 2.0 * m * a * b + nn * b * b)
}
