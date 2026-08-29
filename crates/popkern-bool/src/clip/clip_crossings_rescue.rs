//! Sampled-transition crossing rescue (split out of `clip_uv.rs` for
//! the 500-line ratchet).

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_uv::UvCrossing;

/// Sampled crossing scan: march the query curve densely (adaptive,
/// chords bounded by the trim scale), classify each sample with the
/// face containment predicate, and bisect every in/out transition to
/// a parameter-accurate crossing. Honest by construction — each
/// sample is normalized to the boundary's period branch as a POINT
/// (points cannot tear the way per-control-point curve folding
/// does), and `point_inside_face_uv` is the same predicate the rest
/// of the kernel trusts. Used only as a rescue when the precise pass
/// verifies to nothing.
pub(super) fn sampled_uv_crossings(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    boundary_nurbs: &[(popkern_core::topology::HalfEdgeHandle, NurbsCurve)],
    period: f64,
    tolerance: &Tolerance,
) -> Vec<UvCrossing> {
    use popkern_core::geometry::curve::CurveEval as _;
    let (t0, t1) = curve_uv.parameter_range();
    if !(t1 - t0).is_finite() || t1 - t0 <= 0.0 {
        return Vec::new();
    }
    const NB: usize = 32;
    // ADAPTIVE query polyline: a uniform sampling under-resolves the
    // starvation geometry this rescue exists for — the gusset-web
    // ellipse plunges through the whole 40-unit trim inside one of
    // 512 uniform chords (measured: the rescue found nothing and the
    // phantom window survived). Subdivide until each chord is short
    // relative to the TRIM scale, so no in-face span can hide inside
    // a chord.
    let (b_lu, b_hu, trim_scale) = {
        let (mut lu, mut hu, mut lv, mut hv) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for (_, pc) in boundary_nurbs {
            for p in &pc.control_points {
                lu = lu.min(p.x);
                hu = hu.max(p.x);
                lv = lv.min(p.y);
                hv = hv.max(p.y);
            }
        }
        if lu > hu {
            return Vec::new();
        }
        (
            lu,
            hu,
            ((hu - lu).powi(2) + (hv - lv).powi(2)).sqrt().max(1e-9),
        )
    };
    let max_chord = trim_scale * 0.25;
    let mut q: Vec<(f64, f64, f64)> = Vec::new();
    // Iterative stack subdivision, capped to stay bounded on
    // pathological curves.
    const MAX_POINTS: usize = 65_536;
    type TSample = (f64, f64, f64); // (t, u, v)
    let eval = |t: f64| -> TSample {
        let p = curve_uv.point_at(t);
        (t, p.x, p.y)
    };
    let mut stack: Vec<(TSample, TSample, u32)> = Vec::new();
    // Seed with 64 uniform spans so closed curves cannot collapse.
    let seeds: Vec<TSample> = (0..=64)
        .map(|i| eval(t0 + (t1 - t0) * i as f64 / 64.0))
        .collect();
    for w in seeds.windows(2).rev() {
        stack.push((w[0], w[1], 0));
    }
    q.push(seeds[0]);
    while let Some((a, b, depth)) = stack.pop() {
        let chord = ((b.1 - a.1).powi(2) + (b.2 - a.2).powi(2)).sqrt();
        if chord > max_chord && depth < 16 && q.len() + stack.len() < MAX_POINTS {
            let m = eval(0.5 * (a.0 + b.0));
            stack.push((m, b, depth + 1));
            stack.push((a, m, depth + 1));
        } else {
            q.push(b);
        }
    }
    if popkern_core::trace::trace_enabled("ring-crossing") {
        let (mut lu, mut hu, mut lv, mut hv) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for &(_, u, v) in &q {
            lu = lu.min(u);
            hu = hu.max(u);
            lv = lv.min(v);
            hv = hv.max(v);
        }
        eprintln!(
            "      [xing-scan] {} pts u[{lu:.2},{hu:.2}] v[{lv:.2},{hv:.2}] trim_scale={trim_scale:.2}",
            q.len()
        );
    }
    // Containment per sample; each sample's u is normalized into the
    // boundary's period cell as an independent POINT.
    let b_center_u = 0.5 * (b_lu + b_hu);
    let norm_u = |u: f64| -> f64 {
        if period > 0.0 && period.is_finite() {
            u - period * ((u - b_center_u) / period).round()
        } else {
            u
        }
    };
    let inside_at = |t: f64| -> bool {
        let p = curve_uv.point_at(t);
        crate::split_inner_loop::point_inside_face_uv(
            brep,
            face,
            Point3::new(norm_u(p.x), p.y, 0.0),
            tolerance,
        )
    };
    let flags: Vec<bool> = q.iter().map(|&(t, _, _)| inside_at(t)).collect();
    let mut found: Vec<UvCrossing> = Vec::new();
    for i in 0..q.len().saturating_sub(1) {
        if flags[i] == flags[i + 1] {
            continue;
        }
        // Bisect the transition to a parameter-accurate crossing.
        let (mut ta, mut tb) = (q[i].0, q[i + 1].0);
        let side_a = flags[i];
        for _ in 0..48 {
            let tm = 0.5 * (ta + tb);
            if inside_at(tm) == side_a {
                ta = tm;
            } else {
                tb = tm;
            }
        }
        let t_hit = 0.5 * (ta + tb);
        let ph = curve_uv.point_at(t_hit);
        let (phu, phv) = (norm_u(ph.x), ph.y);
        // Attribute the crossing to the nearest boundary pcurve
        // station (period-normalized u distance).
        let mut best: Option<(f64, popkern_core::topology::HalfEdgeHandle, f64)> = None;
        for (he, pc) in boundary_nurbs {
            let (bp0, bp1) = pc.parameter_range();
            if !(bp1 - bp0).is_finite() || (bp1 - bp0).abs() < 1e-15 {
                continue;
            }
            for i in 0..=NB {
                let tb2 = bp0 + (bp1 - bp0) * i as f64 / NB as f64;
                let p = pc.point_at(tb2);
                let du = norm_u(phu - p.x + b_center_u) - b_center_u;
                let d = du.hypot(phv - p.y);
                if best.map(|(bd, _, _)| d < bd).unwrap_or(true) {
                    best = Some((d, *he, tb2));
                }
            }
        }
        let Some((_, he, t_boundary)) = best else {
            continue;
        };
        if popkern_core::trace::trace_enabled("ring-crossing") {
            eprintln!("      [xing-sampled] he {he:?} t2={t_hit:.4} uv=({phu:.4},{phv:.4})");
        }
        let min_sep = (t1 - t0) * 1e-7;
        if found
            .iter()
            .all(|c: &UvCrossing| (c.t_curve - t_hit).abs() > min_sep)
        {
            found.push(UvCrossing {
                t_curve: t_hit,
                boundary_he: he,
                t_boundary,
            });
        }
    }
    found
}
