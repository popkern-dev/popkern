//! Closed records that ride a face's boundary for their whole length.
//!
//! A plane ∩ wall circle that IS the rim of a disc face (a cone's base
//! rim on a coaxial cylinder wall of the same radius, PK-925's union)
//! lies on the disc's boundary everywhere. The crossing finder then
//! reports coincidence noise as crossings: hits at the record's own
//! start, at the rim's vertex, and a spurious graze one chord short of
//! closure. The crossing hull kept the record short by that chord, and
//! the cylinder wall served split in two (measured: the union rotated
//! about y by -0.4, 0.059 rad of the rim missing, five faces where the
//! unrotated serves four; PK-929). Unrotated, the record starts at the
//! rim's vertex and the noise lands at both ends by luck.
//!
//! Such a record lies in the face's closed region by construction, so
//! the clip keeps it whole. The test is geometric: every sample of the
//! record lies within the boundary band of some boundary pcurve.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::topology::{BRep, FaceHandle};

/// Samples along the record.
const RECORD_SAMPLES: usize = 32;
/// Coarse samples per boundary pcurve before the local refinement.
const BOUNDARY_SAMPLES: usize = 64;

/// The uv pcurves of every loop of `face`.
fn boundary_pcurves(brep: &BRep, face: FaceHandle) -> Option<Vec<Pcurve>> {
    let f = brep.faces.get(face)?;
    let mut out = Vec::new();
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let start = brep.loops.get(l)?.half_edge;
        let mut h = start;
        for _ in 0..=brep.half_edges.len() {
            let he = brep.half_edges.get(h)?;
            out.push(he.pcurve.clone()?);
            h = he.next;
            if h == start {
                break;
            }
        }
    }
    Some(out)
}

/// Uv distance from `q` to `pc`: the nearest coarse sample, refined by
/// a ternary search over its neighbouring parameter span.
fn distance_to(pc: &Pcurve, q: popkern_core::geometry::Point3) -> f64 {
    let (a, b) = pc.parameter_range();
    if !(a.is_finite() && b.is_finite()) {
        return f64::INFINITY;
    }
    let d = |t: f64| {
        let p = pc.point_at(t);
        ((p.x - q.x).powi(2) + (p.y - q.y).powi(2)).sqrt()
    };
    let step = (b - a) / BOUNDARY_SAMPLES as f64;
    let best = (0..=BOUNDARY_SAMPLES)
        .map(|i| a + step * i as f64)
        .min_by(|&s, &t| d(s).total_cmp(&d(t)))
        .unwrap_or(a);
    let (mut lo, mut hi) = ((best - step).max(a), (best + step).min(b));
    for _ in 0..40 {
        let (m1, m2) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if d(m1) < d(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    d(0.5 * (lo + hi))
}

/// Whether every sample of `curve_uv` lies within the clip's boundary
/// band of `face`'s boundary (see the module docs). A face with a
/// missing boundary pcurve never qualifies.
pub(super) fn rides_boundary_everywhere(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
) -> bool {
    let Some(boundary) = boundary_pcurves(brep, face) else {
        return false;
    };
    if boundary.is_empty() {
        return false;
    }
    let band = super::clip_uv::uv_boundary_band(brep, curve_uv);
    let (t0, t1) = curve_uv.parameter_range();
    (0..RECORD_SAMPLES).all(|k| {
        let q = curve_uv.point_at(t0 + (t1 - t0) * (k as f64 + 0.5) / RECORD_SAMPLES as f64);
        boundary.iter().any(|pc| distance_to(pc, q) < band)
    })
}
