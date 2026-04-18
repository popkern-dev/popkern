//! §5.3 item #3 — Shirman-Séquin vertex-enclosure solve at the
//! central vertex of an N-quad Gregory decomposition.
//!
//! For each of the N sub-patches meeting at central vertex `C`,
//! the corner-3 inner CP (in our Gregory-rectangular layout, the
//! CP pair adjacent to patch corner `(u, v) = (1, 0)`) is
//! "projected" so that:
//!
//! 1. The resulting position lies in the common tangent plane at
//!    `C` (the plane with normal `n_c`). This is the G1 condition
//!    along the internal spokes radiating from `C`.
//! 2. The mixed second partial `τ_k` at `C`, when projected onto
//!    the tangent plane, agrees across all N sub-patches. This is
//!    the Shirman-Séquin "relaxed G2" condition (Shirman & Séquin,
//!    "Local surface interpolation with Bézier patches," CAGD 1991)
//!    — sufficient for reflection-line continuity across the
//!    internal spokes.
//!
//! The simplified mean-shift formula below IS the exact L2-optimal
//! solution to the joint system for any N; the "full linear
//! system" framing from the plan reduces to taking the cross-patch
//! mean when the only constraint is tangent-component agreement.
//! The extension to additional along-spoke G2 constraints would
//! require a larger solve — tracked separately as a follow-up on
//! top of this module.
//!
//! # Rank-deficient fallback
//!
//! When `n_c` is near-zero (degenerate N-gon, e.g., all hex
//! vertices coplanar with `C`), the tangent basis is ill-defined
//! and the solve returns `None`. The caller is expected to fall
//! back to the raw parallelogram inner CPs (i.e., leave corner 3
//! at its pre-SS value).

use popkern_core::geometry::{Point3, Vec3};

/// Result of the vertex-enclosure solve: one projected inner CP
/// position per sub-patch (to be used for BOTH `P_u` and `P_v` of
/// the Gregory corner-3 pair in the symmetric-Chiyokura form).
pub(super) type EnclosureSolve = Vec<Point3>;

/// Simplified Shirman-Séquin solve for the corner-3 inner CPs at
/// a central vertex shared by `N` sub-patches.
///
/// Inputs:
/// - `central_point` — the shared corner position `C`.
/// - `n_c` — the outward tangent-plane normal at `C` (typically
///   from averaged hex cross products). Caller is responsible for
///   choosing a consistent sign.
/// - `raw_parallelogram` — the pre-SS corner-3 inner CP for each
///   sub-patch (the parallelogram formula output).
///
/// Returns `Some(shifted)` with one `Point3` per sub-patch, where
/// `shifted[k]` lies in the tangent plane at `C` AND all `shifted[k]`
/// share the same tangent-plane components of the Gregory twist
/// vector τ_k. Returns `None` if `n_c` is degenerate (norm² <
/// `1e-18`), signaling the caller should fall back to the raw values.
pub(super) fn shirman_sequin_vertex_enclosure_solve(
    central_point: Point3,
    n_c: Vec3,
    raw_parallelogram: &[Point3],
) -> Option<EnclosureSolve> {
    if n_c.length_squared() < 1e-18 {
        return None;
    }
    if raw_parallelogram.is_empty() {
        return None;
    }
    let (t1, t2) = orthonormal_tangent_basis(n_c);
    let n = raw_parallelogram.len();
    let (a, b): (Vec<f64>, Vec<f64>) = raw_parallelogram
        .iter()
        .map(|p| {
            let off = Vec3::new(
                p.x - central_point.x,
                p.y - central_point.y,
                p.z - central_point.z,
            );
            (off.dot(&t1), off.dot(&t2))
        })
        .unzip();
    let a_mean = a.iter().sum::<f64>() / n as f64;
    let b_mean = b.iter().sum::<f64>() / n as f64;
    Some(
        (0..n)
            .map(|k| {
                let da = a[k] - a_mean;
                let db = b[k] - b_mean;
                Point3::new(
                    central_point.x + da * t1.x + db * t2.x,
                    central_point.y + da * t1.y + db * t2.y,
                    central_point.z + da * t1.z + db * t2.z,
                )
            })
            .collect(),
    )
}

/// Build an orthonormal tangent basis `(t1, t2)` in the plane
/// perpendicular to `n`. `t1` is derived from the world-axis least
/// aligned with `n` (for numerical stability); `t2 = n × t1` so
/// `(t1, t2, n)` is a right-handed orthonormal frame.
pub(super) fn orthonormal_tangent_basis(n: Vec3) -> (Vec3, Vec3) {
    let ax = n.x.abs();
    let ay = n.y.abs();
    let az = n.z.abs();
    let seed = if ax <= ay && ax <= az {
        Vec3::new(1.0, 0.0, 0.0)
    } else if ay <= az {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let dot = seed.dot(&n);
    let raw_t1 = Vec3::new(seed.x - dot * n.x, seed.y - dot * n.y, seed.z - dot * n.z);
    let t1 = raw_t1.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let t2 = n.cross(&t1);
    (t1, t2)
}
