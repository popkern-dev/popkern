//! Dense-sampled curve projection (split from `clip_endpoint_far.rs`
//! under the 500-line structural budget).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval as _;

/// Dense-sampled projection of `q` onto a record's 3D curve: global
/// best over 256 samples, then a locally-capped Newton polish. The
/// generic 5-start `closest_point` walks the wrong lobe of a long
/// helical record whose distance function has many near-minima
/// (measured on the angled tap: a seam-window piece near z=7.8
/// projected into the exit region at z=8.0 half a curve away, and
/// the far re-trim rebuilt the piece over the wrong stretch).
pub(crate) fn project_onto_curve_dense(
    curve: &popkern_core::geometry::nurbs_curve::NurbsCurve,
    q: &Point3,
) -> Option<f64> {
    let (t0, t1) = curve.parameter_range();
    if !(t1 - t0).is_finite() || (t1 - t0).abs() < 1e-30 {
        return None;
    }
    let n = 256usize;
    let mut bt = t0;
    let mut bd = f64::MAX;
    for k in 0..=n {
        let t = t0 + (t1 - t0) * k as f64 / n as f64;
        let d = curve.point_at(t).distance_to(q);
        if d < bd {
            bd = d;
            bt = t;
        }
    }
    let step_cap = ((t1 - t0) / n as f64).abs();
    let mut t = bt;
    for _ in 0..20 {
        let p = curve.point_at(t);
        let d1 = curve.derivative_at(t, 1);
        let d2 = curve.derivative_at(t, 2);
        let diff = p - *q;
        let f = diff.dot(&d1);
        let fp = d1.dot(&d1) + diff.dot(&d2);
        if fp.abs() < 1e-18 {
            break;
        }
        let dt = (-f / fp).clamp(-step_cap, step_cap);
        t = (t + dt).clamp(t0.min(t1), t0.max(t1));
        if dt.abs() < 1e-13 {
            break;
        }
    }
    Some(t)
}
