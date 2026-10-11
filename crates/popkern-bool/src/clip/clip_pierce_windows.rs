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
use crate::clip_endpoint_refine::trim_boundary_edges;

/// Exact crossing parameters of `curve` with `face_x`'s trim
/// boundary: every boundary edge of the face intersected with the
/// OTHER face's surface, junctions projected onto the record and kept
/// when they lie ON it. Each parameter carries its exact junction
/// POINT (PK-472): a piece end clipped at the parameter alone
/// evaluates on the record's FITTED curve and inherits its lateral
/// slop (~6e-4 measured on the m6 nut's mouth corner — two records
/// meeting at one junction minted nodes that far apart, and the
/// assembled shell failed validation on the gap). `None` when no
/// junction lands on the record (the tuned machinery keeps full
/// authority there).
pub(super) fn pierce_crossing_params(
    brep_x: &BRep,
    face_x: FaceHandle,
    other_brep: &BRep,
    other_face: FaceHandle,
    curve: &CurveOnSurface,
    cache: &mut PierceInventoryCache,
    x_is_a: bool,
) -> Option<Vec<(f64, popkern_core::geometry::Point3)>> {
    let tolerance = &brep_x.tolerance;
    let fit_gate = tolerance.linear() * 2000.0;
    let (t_min, t_max) = curve.curve_3d.parameter_range();
    let span = (t_max - t_min).abs();
    let points = cache.junction_points(brep_x, face_x, other_brep, other_face, x_is_a);
    let mut ts: Vec<(f64, popkern_core::geometry::Point3)> = Vec::new();
    for point in points {
        let Some(t) = project_onto_curve_dense(&curve.curve_3d, &point) else {
            continue;
        };
        if curve.curve_3d.point_at(t).distance_to(&point) > fit_gate {
            continue;
        }
        // The curve ends are already junction-cut; a crossing at
        // an end adds no interval boundary.
        if !curve.is_closed && ((t - t_min).abs() < span * 1e-4 || (t_max - t).abs() < span * 1e-4)
        {
            continue;
        }
        if !ts.iter().any(|&(q, _)| (q - t).abs() < span * 1e-6) {
            ts.push((t, point));
        }
    }
    if ts.is_empty() {
        return None;
    }
    ts.sort_by(|a, b| a.0.total_cmp(&b.0));
    Some(ts)
}

/// Per-face-pair memo of the pierce junction inventory: for each
/// (side, face_x, other_face), the refined edge×surface junction
/// POINTS of `face_x`'s trim boundary against the other face's
/// surface. The inventory is RECORD-INDEPENDENT, but the tuned clip
/// recomputed it per record — a multi-wrap thread flank carries one
/// record per wrap, so the same boundary-coil × shank intersections
/// (each a dense solve over an O(turns)-span coil edge, plus
/// triple-point refinement) were re-derived ~turns times per pair:
/// the intersection stage's ~N^1.8 term (PK-511, stack-sampled:
/// `pierce_crossing_params` → `intersect_surface_curve` under
/// `clip_and_classify`). One cache lives per face-pair task; a memo
/// of a deterministic function, so it cannot change any answer.
#[derive(Default)]
pub(crate) struct PierceInventoryCache {
    points: std::collections::BTreeMap<
        (bool, FaceHandle, FaceHandle),
        Vec<popkern_core::geometry::Point3>,
    >,
}

impl PierceInventoryCache {
    /// The refined junction points of `face_x`'s boundary edges
    /// against `other_face`'s surface (see the struct doc). Computed
    /// once per (side, face pair); the per-record projection and
    /// gating stay in [`pierce_crossing_params`].
    fn junction_points(
        &mut self,
        brep_x: &BRep,
        face_x: FaceHandle,
        other_brep: &BRep,
        other_face: FaceHandle,
        x_is_a: bool,
    ) -> Vec<popkern_core::geometry::Point3> {
        if let Some(points) = self.points.get(&(x_is_a, face_x, other_face)) {
            return points.clone();
        }
        let points = compute_junction_points(brep_x, face_x, other_brep, other_face);
        self.points
            .insert((x_is_a, face_x, other_face), points.clone());
        points
    }
}

/// The pair-level inventory body (formerly inline in
/// [`pierce_crossing_params`], moved verbatim minus the per-record
/// projection): every boundary edge of `face_x` intersected with
/// `other_face`'s surface, each crossing refined to the exact triple
/// point (PK-480 — the 1D solve's position slides along a
/// near-tangent edge by fit-slop/sin(crossing angle); the
/// three-surface Newton removes the amplification, see
/// `clip_triple_junction`).
fn compute_junction_points(
    brep_x: &BRep,
    face_x: FaceHandle,
    other_brep: &BRep,
    other_face: FaceHandle,
) -> Vec<popkern_core::geometry::Point3> {
    let mut points = Vec::new();
    let Some(other_surface) = other_brep
        .faces
        .get(other_face)
        .and_then(|f| f.surface.as_ref())
    else {
        return points;
    };
    let other_surf = &other_surface.geometry;
    let ow = crate::intersect::bounded_parameter_range(other_brep, other_face, other_surf);
    let tolerance = &brep_x.tolerance;
    for eh in trim_boundary_edges(brep_x, face_x) {
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
            points.push(crate::clip_triple_junction::refine_to_triple_point(
                brep_x, eh, other_surf, *point, tolerance,
            ));
        }
    }
    points
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
            trim_boundary_edges(brep, face)
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
        let mut a_yes = 0i32;
        let mut b_yes = 0i32;
        for k in 1..=5 {
            let t = t_lo + (t_hi - t_lo) * k as f64 / 6.0;
            let in_a = region_contains(
                brep_a,
                face_a,
                fold(curve.curve_uv1.point_at(t), &fold_a),
                &fold_a,
            );
            let in_b = region_contains(
                brep_b,
                face_b,
                fold(curve.curve_uv2.point_at(t), &fold_b),
                &fold_b,
            );
            decisive += 1;
            a_yes += i32::from(in_a);
            b_yes += i32::from(in_b);
            if in_a && in_b {
                votes += 1;
            }
        }
        if popkern_core::trace::trace_enabled("clip-far") {
            let mid = curve.curve_3d.point_at(0.5 * (t_lo + t_hi));
            let uv1 = curve.curve_uv1.point_at(0.5 * (t_lo + t_hi));
            let uv2 = curve.curve_uv2.point_at(0.5 * (t_lo + t_hi));
            eprintln!(
                "    [pierce-vote] {face_a:?}x{face_b:?} t[{t_lo:.3},{t_hi:.3}] mid ({:.3},{:.3},{:.3}) uv1 ({:.3},{:.3}) uv2 ({:.3},{:.3}) a {a_yes}/5 b {b_yes}/5 -> {}",
                mid.x,
                mid.y,
                mid.z,
                uv1.x,
                uv1.y,
                uv2.x,
                uv2.y,
                votes * 2 > decisive
            );
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

/// A vote sample's face-region containment. On a chart with no period
/// the decision is EXACT near curved rims ([`point_in_face_region_exact`]):
/// the tuned test polygonises each boundary pcurve at 8 stations, and a
/// boss cap of radius 42.5 became an octagon whose chords cut 3.2 inside
/// its rim. Votes 0.18 inside the rim (a ring-wall record between the
/// cap's top and the web notch, the boss cap 0.1 inside the ring's
/// extreme) read Outside, the override dropped that window, and the
/// union stranded the cap's top sliver (PK-780). Periodic charts and
/// rings without pcurves keep the tuned test.
///
/// [`point_in_face_region_exact`]: crate::classify_uv_loops::point_in_face_region_exact
fn region_contains(
    brep: &BRep,
    face: FaceHandle,
    uv: popkern_core::geometry::Point3,
    chart: &(Option<f64>, Option<f64>, (f64, f64)),
) -> bool {
    if chart.0.is_none()
        && chart.1.is_none()
        && let Some(inside) =
            crate::classify_uv_loops::point_in_face_region_exact(brep, face, uv.x, uv.y)
    {
        return inside;
    }
    crate::split_inner_loop::point_inside_face_region_uv(brep, face, uv, &brep.tolerance)
}

/// The PIERCE-AUTHORITY override applied to a tuned joint window set
/// (extracted from `clip_windows.rs` under the 500-line budget):
/// near-tangent boundary crossings defeat the pcurve crossing
/// machinery — the exact edge×surface junction inventory re-derives
/// the windows and overrides the tuned result only where the two
/// disagree. Scoped by the caller to transversal OPEN records with no
/// wrapped tuned window. Returns the (possibly overridden) joint
/// windows and the certified junction points backing pierce-derived
/// window bounds (PK-472: piece ends clipped at those bounds pin to
/// the junction, not the fitted curve's lateral slop).
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(super) fn apply_pierce_authority(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
    t_min: f64,
    t_max: f64,
    joint: Vec<(f64, f64)>,
    trace: bool,
    cache: &mut PierceInventoryCache,
) -> (Vec<(f64, f64)>, Vec<(f64, popkern_core::geometry::Point3)>) {
    let mut joint = joint;
    let mut junction_points: Vec<(f64, popkern_core::geometry::Point3)> = Vec::new();
    // PIERCE AUTHORITY (see `clip_pierce_windows`): near-tangent
    // boundary crossings defeat the pcurve crossing machinery — the
    // exact edge×surface junction inventory re-derives the windows,
    // and overrides the tuned result only where the two disagree.
    // Scoped to OPEN records with no wrapped tuned window: the
    // wrapped/closed bookkeeping keeps its tuned path.
    // Transversal records only — a flush/coincident trace RUNS ALONG
    // boundaries, its pierce inventory is everywhere on it, and the
    // interval containment votes sit inside boundary slop (measured:
    // the notched-cube saddle fixture's flush records got mangled
    // windows and the subtraction went non-manifold).
    {
        let span = (t_max - t_min).abs();
        let ts_a = pierce_crossing_params(brep_a, face_a, brep_b, face_b, curve, cache, true);
        let ts_b = pierce_crossing_params(brep_b, face_b, brep_a, face_a, curve, cache, false);
        if ts_a.is_some() || ts_b.is_some() {
            let mut junctions = ts_a.unwrap_or_default();
            junctions.extend(ts_b.unwrap_or_default());
            junctions.sort_by(|a, b| a.0.total_cmp(&b.0));
            junctions.dedup_by(|b, a| (b.0 - a.0).abs() < span * 1e-6);
            let mut ts: Vec<f64> = junctions.iter().map(|&(t, _)| t).collect();
            if popkern_core::trace::trace_enabled("clip-far") {
                let (t0m, t1m) = curve.curve_3d.parameter_range();
                let p0 = curve.curve_3d.point_at(t0m);
                let p1 = curve.curve_3d.point_at(t1m);
                eprintln!(
                    "    [pierce-ts] record ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}): ts {:?}",
                    p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, ts
                );
            }
            let cleared = intervals_clear_of_boundary(brep_a, face_a, brep_b, face_b, curve, &ts);
            if !cleared {
                ts.clear();
            }
            let pw = pierce_windows(brep_a, face_a, brep_b, face_b, curve, &ts);
            let close = |x: f64, y: f64| (x - y).abs() <= span * 1e-3;
            let differs = pw.len() != joint.len()
                || pw
                    .iter()
                    .zip(joint.iter())
                    .any(|(x, y)| !close(x.0, y.0) || !close(x.1, y.1));
            // A pierce result that deletes EVERYTHING is never
            // trusted: the measured defect class needs windows
            // inverted or split, not erased, and an all-outside
            // verdict is the signature of a containment context the
            // voter cannot read (extended-base pcurve branches sit
            // whole chart-widths from the boundary bbox — measured on
            // the bore-mouth family, where a correct tuned window was
            // emptied).
            // Pins arm whenever CERTIFIED junctions exist — matching
            // stays EXACT (span x 1e-6), so only window bounds that
            // came from the junction machinery ever pin. Two face
            // pairs' chains solving one triple junction then land on
            // the SAME exact point instead of two nodes a lateral
            // slop apart (measured, PK-465: the angled tap's mouth
            // minted 5.4e-4 record-record duplicates the weld tiers
            // cannot own — 7 stranded sides). Slop-gated arming was
            // measured harmful (re-strands) and stays out.
            if cleared && !ts.is_empty() {
                junction_points = junctions.clone();
            }
            if cleared && !ts.is_empty() && !pw.is_empty() && differs {
                if trace {
                    let show = |ws: &[(f64, f64)]| {
                        ws.iter()
                            .map(|(e, x)| format!("({e:.3},{x:.3})"))
                            .collect::<Vec<_>>()
                            .join(",")
                    };
                    eprintln!(
                        "  [uvclip] pierce-authority override: tuned [{}] -> exact [{}] from {} crossing(s)",
                        show(&joint),
                        show(&pw),
                        ts.len()
                    );
                }
                joint = pw;
            }
        }
    }
    (joint, junction_points)
}
