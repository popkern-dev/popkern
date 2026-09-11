//! Window-set algebra and record transversality for the clip
//! (split from `clip_windows.rs` under the 500-line file budget).

use popkern_core::topology::{BRep, FaceHandle};

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;

/// Measured transversality of an intersection record: the two faces'
/// unit normals must disagree by at least ~3 degrees at every probed
/// record sample. Coincident and tangential pairs (the flush
/// boolean's coplanar overlap class) fail this and keep the tuned
/// crossing-hull clip.
pub(super) fn record_is_transversal(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    let (Some(surface_a), Some(surface_b)) = (
        brep_a.faces.get(face_a).and_then(|f| f.surface.as_ref()),
        brep_b.faces.get(face_b).and_then(|f| f.surface.as_ref()),
    ) else {
        return false;
    };
    let (a_lo, a_hi) = curve.curve_uv1.parameter_range();
    let (b_lo, b_hi) = curve.curve_uv2.parameter_range();
    const MIN_NORMAL_CROSS: f64 = 5.0e-2;
    for fraction in [0.25, 0.5, 0.75] {
        let uv_a = curve.curve_uv1.point_at(a_lo + fraction * (a_hi - a_lo));
        let uv_b = curve.curve_uv2.point_at(b_lo + fraction * (b_hi - b_lo));
        let normal_a = surface_a.geometry.normal_at(uv_a.x, uv_a.y);
        let normal_b = surface_b.geometry.normal_at(uv_b.x, uv_b.y);
        let (Some(unit_a), Some(unit_b)) = (normal_a.normalized(), normal_b.normalized()) else {
            return false;
        };
        let cross = unit_a.cross(&unit_b).length();
        if cross < MIN_NORMAL_CROSS {
            if popkern_core::trace::trace_enabled("clip-far") {
                use popkern_core::geometry::curve::CurveEval as _;
                let (t0, t1) = curve.curve_3d.parameter_range();
                let pm = curve.curve_3d.point_at(t0 + fraction * (t1 - t0));
                eprintln!(
                    "    [transversal] record ({:.3},{:.3},{:.3})-sample f{fraction} cross {cross:.3e} -> NOT transversal",
                    pm.x, pm.y, pm.z
                );
            }
            return false;
        }
    }
    true
}

/// Intersect two in-face window sets on the record's parameter line —
/// on its parameter CIRCLE for closed records, where either side's
/// windows may run through the seam (exit beyond the range end) and a
/// single window pair can overlap in two separate arcs.
pub(super) fn intersect_window_sets(
    wa: &[(f64, f64)],
    wb: &[(f64, f64)],
    t_min: f64,
    t_max: f64,
    is_closed: bool,
) -> Vec<(f64, f64)> {
    let period = t_max - t_min;
    let full = |w: &(f64, f64)| w.1 - w.0 >= period - 1e-12;
    // A full-range window is no constraint at all: the other set
    // passes through unchanged (in particular a seam-running window
    // survives as ONE wrapped window, not two split arcs).
    if let [w] = wa
        && full(w)
    {
        return wb.to_vec();
    }
    if let [w] = wb
        && full(w)
    {
        return wa.to_vec();
    }
    let shifts: &[f64] = if is_closed { &[-1.0, 0.0, 1.0] } else { &[0.0] };
    let mut out: Vec<(f64, f64)> = Vec::new();
    for a in wa {
        for b in wb {
            for &s in shifts {
                let e = a.0.max(b.0 + s * period);
                let x = a.1.min(b.1 + s * period);
                if x - e <= 1e-15 {
                    continue;
                }
                let (mut e, mut x) = (e, x);
                if is_closed {
                    let k = ((e - t_min) / period).floor();
                    e -= k * period;
                    x -= k * period;
                }
                out.push((e, x));
            }
        }
    }
    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let eps = period.abs().max(1e-12) * 1e-9;
    out.dedup_by(|b, a| (b.0 - a.0).abs() < eps && (b.1 - a.1).abs() < eps);
    out
}

/// Subtract window set `wb` from `wa` on the record's parameter
/// line: the spans inside `wa` but outside every `wb` window. OPEN
/// records only (residue extraction is scoped to transversal open
/// records — see `clip_curve_uv_space`).
pub(super) fn subtract_window_sets(wa: &[(f64, f64)], wb: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = wa.to_vec();
    for b in wb {
        let mut next = Vec::with_capacity(out.len() + 1);
        for a in &out {
            if b.1 <= a.0 || b.0 >= a.1 {
                next.push(*a);
                continue;
            }
            if b.0 > a.0 {
                next.push((a.0, b.0));
            }
            if b.1 < a.1 {
                next.push((b.1, a.1));
            }
        }
        out = next;
    }
    out
}

/// the clip trace topic preamble: the record's 3D/UV extents and the B
/// face's boundary bbox.
pub(super) fn trace_record_geometry(
    brep_b: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
) {
    let (a3, b3) = curve.curve_3d.parameter_range();
    let (a1, b1) = curve.curve_uv1.parameter_range();
    let p3a = curve.curve_3d.point_at(a3);
    let p3b = curve.curve_3d.point_at(b3);
    let p1a = curve.curve_uv1.point_at(a1);
    let p1b = curve.curve_uv1.point_at(b1);
    let (a2, b2) = curve.curve_uv2.parameter_range();
    let p2a = curve.curve_uv2.point_at(a2);
    let p2b = curve.curve_uv2.point_at(b2);
    eprintln!(
        "    [uvspace] {face_a:?}x{face_b:?} 3d t[{a3:.1},{b3:.1}] ({:.1},{:.1},{:.1})->({:.1},{:.1},{:.1}) | uv1 t[{a1:.1},{b1:.1}] ({:.1},{:.1})->({:.1},{:.1}) | uv2 ({:.3},{:.1})->({:.3},{:.1})",
        p3a.x,
        p3a.y,
        p3a.z,
        p3b.x,
        p3b.y,
        p3b.z,
        p1a.x,
        p1a.y,
        p1b.x,
        p1b.y,
        p2a.x,
        p2a.y,
        p2b.x,
        p2b.y
    );
    let ring = crate::clip::face_uv_boundary(brep_b, face_b);
    let mut lo = (f64::MAX, f64::MAX);
    let mut hi = (f64::MIN, f64::MIN);
    for (_, pc) in ring.iter().flatten() {
        let (a, b) = pc.parameter_range();
        for i in 0..4 {
            let q = pc.point_at(a + (b - a) * i as f64 / 3.0);
            lo = (lo.0.min(q.x), lo.1.min(q.y));
            hi = (hi.0.max(q.x), hi.1.max(q.y));
        }
    }
    eprintln!(
        "    [uvspace] B-poly bbox ({:.3},{:.1})..({:.3},{:.1})",
        lo.0, lo.1, hi.0, hi.1
    );
}
