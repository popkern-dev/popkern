//! Pierce-authoritative clip windows.
//!
//! A record crossing a face boundary NEAR-TANGENT defeats the pcurve
//! crossing machinery: the crossing parameters displace by
//! 1/sin(crossing angle), crossings get missed or mis-ordered, and a
//! window verdict can INVERT (measured, PK-465: on the angled tap's
//! tilted bore wall, a helical wall×flank turn weaves across the
//! mouth-rim sinusoid — the tuned clip kept exactly the above-rim
//! bulge and dropped both in-face stretches, and every exit arc
//! became a half-claimed strand).
//!
//! The record's true boundary crossings are available EXACTLY: every
//! crossing with face X's trim lies on some boundary EDGE of X and on
//! the OTHER face's surface (the record lives on that surface), so
//! the edge×surface junction inventory — 1D solves that stay robust
//! at any crossing angle — is a complete crossing list. Windows
//! derive from it directly: sort the junction parameters along the
//! record, then vote each interval in/out by interior region
//! containment samples (well away from the boundary, where
//! containment is decisive).
//!
//! Authority is scoped to disagreement: when the tuned machinery's
//! windows already match the pierce-derived ones, the tuned result is
//! kept byte-identical.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip_endpoint_far::project_onto_curve_dense;
use crate::clip_endpoint_refine::boundary_edges;

/// Exact crossing parameters of `curve` with `face_x`'s trim
/// boundary: every boundary edge of the face intersected with the
/// OTHER face's surface, junctions projected onto the record and kept
/// when they lie ON it. `None` when no junction lands on the record
/// (the tuned machinery keeps full authority there).
pub(super) fn pierce_crossing_params(
    brep_x: &BRep,
    face_x: FaceHandle,
    other_brep: &BRep,
    other_face: FaceHandle,
    curve: &CurveOnSurface,
) -> Option<Vec<f64>> {
    let other_surf = &other_brep.faces.get(other_face)?.surface.as_ref()?.geometry;
    let ow = crate::intersect::bounded_parameter_range(other_brep, other_face, other_surf);
    let tolerance = &brep_x.tolerance;
    let fit_gate = tolerance.linear() * 2000.0;
    let (t_min, t_max) = curve.curve_3d.parameter_range();
    let span = (t_max - t_min).abs();
    let mut ts: Vec<f64> = Vec::new();
    for eh in boundary_edges(brep_x, face_x) {
        let Some(edge) = brep_x.edges.get(eh) else {
            continue;
        };
        let Some(c) = edge.curve.as_ref() else {
            continue;
        };
        let range = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
        let geometry: CurveGeometry = c.geometry.clone();
        let Ok(result) = popkern_core::geometry::intersection::intersect_surface_curve(
            other_surf, ow.0, ow.1, &geometry, range, tolerance,
        ) else {
            continue;
        };
        if popkern_core::trace::trace_enabled("clip-far") {
            eprintln!(
                "    [pierce-inv] {face_x:?} edge {eh:?} x other-surf: {} isect(s)",
                result.intersections.len()
            );
        }
        for isect in &result.intersections {
            let popkern_core::geometry::intersection::SurfaceCurveIntersection::Point {
                point, ..
            } = isect
            else {
                continue;
            };
            let Some(t) = project_onto_curve_dense(&curve.curve_3d, point) else {
                continue;
            };
            if curve.curve_3d.point_at(t).distance_to(point) > fit_gate {
                continue;
            }
            // The curve ends are already junction-cut; a crossing at
            // an end adds no interval boundary.
            if !curve.is_closed
                && ((t - t_min).abs() < span * 1e-4 || (t_max - t).abs() < span * 1e-4)
            {
                continue;
            }
            if !ts.iter().any(|&q| (q - t).abs() < span * 1e-6) {
                ts.push(t);
            }
        }
    }
    if ts.is_empty() {
        return None;
    }
    ts.sort_by(f64::total_cmp);
    Some(ts)
}

/// Clearance certificate arming the pierce override: every interval
/// midpoint must sit CLEAR of both faces' boundary edges in 3D. A
/// record that RUNS ALONG a boundary (a straight bore's rim circle is
/// itself the wall∩top junction curve) has its junction inventory
/// dense on it and every containment sample lands within boundary
/// slop — the votes are coin flips and the override destroyed
/// legitimate windows (measured: the bore-mouth corpus family). A
/// genuinely CROSSING record's interval interiors stand off the
/// boundary by real distance.
pub(super) fn intervals_clear_of_boundary(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
    ts: &[f64],
) -> bool {
    let (t_min, t_max) = curve.curve_3d.parameter_range();
    let clearance = brep_a.tolerance.linear().max(brep_b.tolerance.linear()) * 1000.0;
    let mut bounds: Vec<f64> = Vec::with_capacity(ts.len() + 2);
    bounds.push(t_min);
    bounds.extend(ts.iter().copied());
    bounds.push(t_max);
    let edge_curves: Vec<_> = [(brep_a, face_a), (brep_b, face_b)]
        .into_iter()
        .flat_map(|(brep, face)| {
            boundary_edges(brep, face)
                .into_iter()
                .filter_map(move |eh| brep.edges.get(eh)?.curve.clone())
        })
        .collect();
    for w in bounds.windows(2) {
        if (w[1] - w[0]).abs() < (t_max - t_min).abs() * 1e-6 {
            continue;
        }
        let p = curve.curve_3d.point_at(0.5 * (w[0] + w[1]));
        for c in &edge_curves {
            let Some(t) = (&c.geometry as &dyn CurveEval).closest_point(&p) else {
                continue;
            };
            let tc = t.clamp(c.t_start.min(c.t_end), c.t_start.max(c.t_end));
            if c.geometry.point_at(tc).distance_to(&p) < clearance {
                return false;
            }
        }
    }
    true
}

/// Windows over `curve` from the exact crossing inventory: intervals
/// between consecutive crossings (plus the curve ends for open
/// records), each voted in/out by interior region-containment samples
/// on BOTH faces. Adjacent in-face intervals merge.
pub(super) fn pierce_windows(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
    ts: &[f64],
) -> Vec<(f64, f64)> {
    let (t_min, t_max) = curve.curve_3d.parameter_range();
    let mut bounds: Vec<f64> = Vec::with_capacity(ts.len() + 2);
    bounds.push(t_min);
    bounds.extend(ts.iter().copied());
    bounds.push(t_max);
    // PERIOD FOLD before containment: on a periodic chart the
    // record's pcurve can live on another period branch of the
    // extended base (measured: a bore-wall record carried v = 71
    // against a boundary spanning v ∈ [0,60] and every vote read
    // outside, emptying a correct window). Fold each sample by the
    // chart period toward the face's own boundary bbox mid.
    let fold_for = |brep: &BRep, face: FaceHandle| -> (Option<f64>, Option<f64>, (f64, f64)) {
        let periods = brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| popkern_core::geometry::chart_periods::geometric_chart_periods(&s.geometry))
            .unwrap_or((None, None));
        let ring = crate::clip::face_uv_boundary(brep, face);
        let mut lo = (f64::MAX, f64::MAX);
        let mut hi = (f64::MIN, f64::MIN);
        for (_, pc) in ring.iter().flatten() {
            let (a, b) = pc.parameter_range();
            for i in 0..5 {
                let q = pc.point_at(a + (b - a) * i as f64 / 4.0);
                lo = (lo.0.min(q.x), lo.1.min(q.y));
                hi = (hi.0.max(q.x), hi.1.max(q.y));
            }
        }
        let mid = (0.5 * (lo.0 + hi.0), 0.5 * (lo.1 + hi.1));
        (periods.0, periods.1, mid)
    };
    let fold = |q: popkern_core::geometry::Point3,
                f: &(Option<f64>, Option<f64>, (f64, f64))|
     -> popkern_core::geometry::Point3 {
        let fold1 = |x: f64, period: Option<f64>, mid: f64| -> f64 {
            match period {
                Some(p) if p.is_finite() && p > 0.0 => x - ((x - mid) / p).round() * p,
                _ => x,
            }
        };
        popkern_core::geometry::Point3::new(fold1(q.x, f.0, f.2.0), fold1(q.y, f.1, f.2.1), 0.0)
    };
    let fold_a = fold_for(brep_a, face_a);
    let fold_b = fold_for(brep_b, face_b);
    let inside = |t_lo: f64, t_hi: f64| -> bool {
        let mut votes = 0i32;
        let mut decisive = 0i32;
        for k in 1..=5 {
            let t = t_lo + (t_hi - t_lo) * k as f64 / 6.0;
            let in_a = crate::split_inner_loop::point_inside_face_region_uv(
                brep_a,
                face_a,
                fold(curve.curve_uv1.point_at(t), &fold_a),
                &brep_a.tolerance,
            );
            let in_b = crate::split_inner_loop::point_inside_face_region_uv(
                brep_b,
                face_b,
                fold(curve.curve_uv2.point_at(t), &fold_b),
                &brep_b.tolerance,
            );
            decisive += 1;
            if in_a && in_b {
                votes += 1;
            }
        }
        votes * 2 > decisive
    };
    let mut out: Vec<(f64, f64)> = Vec::new();
    for w in bounds.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if hi - lo < (t_max - t_min).abs() * 1e-6 {
            continue;
        }
        if inside(lo, hi) {
            match out.last_mut() {
                Some(last) if (last.1 - lo).abs() < (t_max - t_min).abs() * 1e-9 => last.1 = hi,
                _ => out.push((lo, hi)),
            }
        }
    }
    out
}
