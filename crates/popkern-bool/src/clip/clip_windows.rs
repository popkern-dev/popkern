//! Window-SET clipping of intersection records to face boundaries.
//!
//! A record curve can enter and leave a face's trimmed domain several
//! times — an analytic sphere×plane circle crosses a cube face's
//! bounds four times, a line through a reflex outline twice per bite.
//! The per-face verdict is therefore a set of in-face parameter
//! windows ([`UvClipResult::Windows`]), and the joint verdict is the
//! pairwise intersection of the two faces' window sets (mod the
//! parameter period for closed records). Every surviving window
//! becomes its own clipped `CurveOnSurface` — collapsing them to one
//! span (the pre-PK-111 behavior) manufactured records that ran
//! THROUGH regions outside the face and over foreign junctions.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_uv::{UvClipResult, clip_curve_to_crossings, clip_curve_wrapped, uv_clip_one_face};

/// Measured transversality of an intersection record: the two faces'
/// unit normals must disagree by at least ~3 degrees at every probed
/// record sample. Coincident and tangential pairs (the flush
/// boolean's coplanar overlap class) fail this and keep the tuned
/// crossing-hull clip.
fn record_is_transversal(
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
        if unit_a.cross(&unit_b).length() < MIN_NORMAL_CROSS {
            return false;
        }
    }
    true
}

/// Intersect two in-face window sets on the record's parameter line —
/// on its parameter CIRCLE for closed records, where either side's
/// windows may run through the seam (exit beyond the range end) and a
/// single window pair can overlap in two separate arcs.
fn intersect_window_sets(
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

/// Clip an intersection record to both faces' trimmed UV domains.
/// Returns EVERY in-face span as its own `CurveOnSurface` (empty =
/// the record misses at least one face entirely). A record judged
/// in-face over its whole closed range survives as a single closed
/// clone.
pub(crate) fn clip_curve_uv_space(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
) -> Vec<CurveOnSurface> {
    let trace = popkern_core::trace::trace_enabled("clip");
    if trace {
        trace_record_geometry(brep_b, face_a, face_b, curve);
    }
    // Open records may take decisive multi-windows only when the
    // record is TRANSVERSAL — the two faces' surface normals disagree
    // along it. Flush/coincident boundary-running traces (parallel
    // normals) keep the tuned crossing-hull contract (PK-114).
    let open_windows_allowed =
        !curve.is_closed && record_is_transversal(brep_a, face_a, brep_b, face_b, curve);
    // Chart-cell route (PK-304): doubly-periodic faces take
    // folded-containment windows (see `chart_cell`); everything else
    // keeps the tuned crossing machinery.
    let cell_clip = super::chart_cell::cell_clip;
    let clip_a = cell_clip(brep_a, face_a, &curve.curve_uv1).unwrap_or_else(|| {
        uv_clip_one_face(
            brep_a,
            face_a,
            &curve.curve_uv1,
            curve.is_closed,
            open_windows_allowed,
        )
    });
    let clip_b = cell_clip(brep_b, face_b, &curve.curve_uv2).unwrap_or_else(|| {
        uv_clip_one_face(
            brep_b,
            face_b,
            &curve.curve_uv2,
            curve.is_closed,
            open_windows_allowed,
        )
    });

    let name = |r: &UvClipResult| match r {
        UvClipResult::NoPcurves => "NoPcurves".to_string(),
        UvClipResult::Outside => "Outside".to_string(),
        UvClipResult::Windows(ws) => format!(
            "Windows[{}]",
            ws.iter()
                .map(|(e, x)| format!("({e:.3},{x:.3})"))
                .collect::<Vec<_>>()
                .join(",")
        ),
    };
    if matches!(clip_a, UvClipResult::Outside) || matches!(clip_b, UvClipResult::Outside) {
        if trace {
            eprintln!(
                "  [uvclip] dropped: A={} B={}",
                name(&clip_a),
                name(&clip_b)
            );
        }
        return Vec::new();
    }
    if matches!(clip_a, UvClipResult::NoPcurves) && matches!(clip_b, UvClipResult::NoPcurves) {
        return vec![curve.clone()];
    }
    if trace {
        eprintln!(
            "  [uvclip] {face_a:?}x{face_b:?}: A={} B={}",
            name(&clip_a),
            name(&clip_b)
        );
    }

    let (t_min, t_max) = curve.curve_3d.parameter_range();
    let windows_of = |r: UvClipResult| -> Vec<(f64, f64)> {
        match r {
            UvClipResult::Windows(ws) => ws,
            // NoPcurves on one side: that face imposes no constraint.
            _ => vec![(t_min, t_max)],
        }
    };
    let wa = windows_of(clip_a);
    let wb = windows_of(clip_b);
    let joint = intersect_window_sets(&wa, &wb, t_min, t_max, curve.is_closed);
    if trace && joint.len() > 1 {
        eprintln!(
            "  [uvclip] {} in-face windows survive the joint clip",
            joint.len()
        );
    }

    let mut out = Vec::new();
    for (t_enter, t_exit) in joint {
        if t_exit - t_enter < 1e-15 {
            continue;
        }
        if curve.is_closed && (t_enter - t_min).abs() < 1e-12 && (t_exit - t_max).abs() < 1e-12 {
            out.push(curve.clone());
            continue;
        }
        let wrapped = t_exit > t_max + 1e-12;
        let extract = |c: &NurbsCurve| -> Option<NurbsCurve> {
            if wrapped {
                clip_curve_wrapped(c, t_enter, t_exit)
            } else {
                clip_curve_to_crossings(c, t_enter, t_exit)
            }
        };
        let (Some(clipped_3d), Some(clipped_uv1), Some(clipped_uv2)) = (
            extract(&curve.curve_3d),
            extract(&curve.curve_uv1),
            extract(&curve.curve_uv2),
        ) else {
            continue;
        };
        out.push(CurveOnSurface {
            curve_3d: clipped_3d,
            curve_uv1: clipped_uv1,
            curve_uv2: clipped_uv2,
            is_closed: false,
        });
    }
    out
}

/// the clip trace topic preamble: the record's 3D/UV extents and the B
/// face's boundary bbox.
fn trace_record_geometry(
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
