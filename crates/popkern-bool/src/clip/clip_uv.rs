use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::intersection::intersect_curves;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::face_uv_boundary;

const TWO_PI: f64 = std::f64::consts::TAU;

pub(super) fn is_near_seam(curve: &NurbsCurve) -> bool {
    let quarter = std::f64::consts::FRAC_PI_4;
    let seven_quarter = 7.0 * quarter;
    curve
        .control_points
        .iter()
        .any(|p| p.x < quarter || p.x > seven_quarter)
}

pub(super) fn shift_curve_u(curve: &NurbsCurve, shift: f64) -> NurbsCurve {
    let shifted_points = curve
        .control_points
        .iter()
        .map(|p| {
            let mut u = p.x + shift;
            u = u.rem_euclid(TWO_PI);
            Point3::new(u, p.y, p.z)
        })
        .collect();
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: shifted_points,
        weights: curve.weights.clone(),
    }
}

/// A crossing between a UV intersection curve and a face's UV boundary.
pub(crate) struct UvCrossing {
    pub t_curve: f64,
    pub boundary_he: popkern_core::topology::HalfEdgeHandle,
    pub t_boundary: f64,
}

/// Sampled crossing scan: march the query curve densely (adaptive,
/// chords bounded by the trim scale), classify each sample with the
/// face containment predicate, and bisect every in/out transition to
/// a parameter-accurate crossing. Honest by construction — each
/// sample is normalized to the boundary's period branch as a POINT
/// (points cannot tear the way per-control-point curve folding
/// does), and `point_inside_face_uv` is the same predicate the rest
/// of the kernel trusts. Used only as a rescue when the precise pass
/// verifies to nothing.
fn sampled_uv_crossings(
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
    if std::env::var_os("PK_XING_TRACE").is_some() {
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
        if std::env::var_os("PK_XING_TRACE").is_some() {
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

/// Find where a UV intersection curve crosses the face's UV boundary.
pub(crate) fn find_uv_crossings(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    tolerance: &Tolerance,
) -> Option<Vec<UvCrossing>> {
    let boundary = face_uv_boundary(brep, face)?;
    let mut crossings = Vec::new();

    let curve_range = curve_uv.parameter_range();
    if (curve_range.1 - curve_range.0).abs() < tolerance.numerical_floor() {
        return Some(crossings);
    }

    let in_periodic_range = |curve: &NurbsCurve| -> bool {
        curve
            .control_points
            .iter()
            .any(|p| p.x > std::f64::consts::PI)
    };
    let boundary_nurbs: Vec<(popkern_core::topology::HalfEdgeHandle, NurbsCurve)> =
        boundary.iter().map(|(h, pc)| (*h, pc.to_nurbs())).collect();
    // The seam-shift dance is ONLY meaningful on periodic surfaces:
    // its "near seam" / "in periodic range" tests are raw u-value
    // thresholds (π/4, π) and fire spuriously on PLANAR faces whose
    // UV coordinates are plain lengths (measured: a 20mm-wide tool
    // wall got its intersection segment shifted by π and clipped to
    // ±5.78 instead of ±10, so the pocket never met the wall
    // corners and the boolean output carried degenerate stubs).
    let surface_is_periodic = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| {
            matches!(
                s.geometry,
                popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }
                    | popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. }
                    | popkern_core::geometry::surface::SurfaceGeometry::Torus { .. }
                    | popkern_core::geometry::surface::SurfaceGeometry::Cone { .. }
            )
        })
        .unwrap_or(false);
    // PERIOD-BRANCH NORMALIZATION: the query curve may live on any
    // whole-period u-branch of a periodic chart (a plane×cylinder
    // marcher parameterized one gusset-web ellipse at u ∈ [2π, 4π]),
    // while the face's pcurves live on the principal branch. UV
    // curve×curve intersection then finds no true crossings — and
    // the ±π seam dance below cannot bridge a full period (measured:
    // the arc-web×underside record clipped to the EXTENDED
    // cylinder's phantom re-crossing instead of the true weld line,
    // and the fanned arm's gusset union went non-manifold). Shift
    // the query curve by whole periods so its u-centre lands nearest
    // the boundary's u-centre; crossing t-parameters are unaffected.
    let curve_uv = &{
        let mut normalized = curve_uv.clone();
        if surface_is_periodic && !boundary_nurbs.is_empty() {
            let period = {
                use popkern_core::geometry::SurfaceEval as _;
                brep.faces
                    .get(face)
                    .and_then(|f| f.surface.as_ref())
                    .and_then(|s| s.geometry.parameter_periodicity().0)
                    .unwrap_or(std::f64::consts::TAU)
            };
            let u_mean = |c: &NurbsCurve| {
                c.control_points.iter().map(|p| p.x).sum::<f64>() / c.control_points.len() as f64
            };
            let boundary_center = boundary_nurbs.iter().map(|(_, pc)| u_mean(pc)).sum::<f64>()
                / boundary_nurbs.len() as f64;
            let k = ((boundary_center - u_mean(&normalized)) / period).round();
            if k != 0.0 {
                normalized = shift_curve_u(&normalized, k * period);
                if std::env::var_os("PK_CLIP_TRACE").is_some() {
                    eprintln!(
                        "    [uvclip-branch] {face:?} query curve shifted by {k:+.0} period(s)"
                    );
                }
            }
        }
        normalized
    };
    let needs_shift = surface_is_periodic
        && (is_near_seam(curve_uv) || boundary_nurbs.iter().any(|(_, pc)| is_near_seam(pc)))
        && (in_periodic_range(curve_uv)
            || boundary_nurbs.iter().any(|(_, pc)| in_periodic_range(pc)));
    let shift = std::f64::consts::PI;

    let collect_crossings = |shifted: bool| -> Vec<UvCrossing> {
        let working_curve = if shifted {
            shift_curve_u(curve_uv, shift)
        } else {
            curve_uv.clone()
        };
        let mut found = Vec::new();
        for (he_handle, boundary_pcurve) in &boundary_nurbs {
            let boundary_range = boundary_pcurve.parameter_range();
            if (boundary_range.1 - boundary_range.0).abs() < tolerance.numerical_floor() {
                continue;
            }

            let twin_he = brep.half_edges[*he_handle].twin;
            if brep.half_edges.contains_key(twin_he)
                && brep.half_edges[twin_he].loop_ == brep.half_edges[*he_handle].loop_
            {
                continue;
            }

            let working_boundary = if shifted {
                shift_curve_u(boundary_pcurve, shift)
            } else {
                boundary_pcurve.clone()
            };

            let result = intersect_curves(
                &CurveGeometry::Nurbs(working_boundary),
                boundary_range,
                &CurveGeometry::Nurbs(working_curve.clone()),
                curve_range,
                tolerance,
            );

            if let Ok(result) = result {
                use popkern_core::geometry::intersection::CurveIntersection;
                for hit in &result.intersections {
                    match hit {
                        CurveIntersection::Point { t1, t2, .. } => {
                            // TEAR CONSISTENCY: the fold wraps each
                            // control point by rem_euclid, which
                            // TEARS a curve spanning the wrap — its
                            // NURBS evaluation then sweeps spurious
                            // bridge paths, and a crossing found on a
                            // bridge has a t-parameter pointing at
                            // unrelated geometry (measured: the
                            // gusset-web ellipse's dance crossings
                            // landed on the extended cylinder's
                            // phantom re-crossing half a period from
                            // the boundary they claimed to cross). A
                            // fold-found crossing is honest only
                            // where the folded QUERY evaluation
                            // equals the per-point fold of the
                            // original evaluation. The check is
                            // query-only: BOUNDARY pcurves may
                            // legitimately encode a seam-spanning
                            // edge "the long way" in the chart, and
                            // tearing them is precisely the healing
                            // the dance exists for; query records
                            // come from the marcher continuously
                            // unwrapped, so a torn query is always a
                            // bridge artifact.
                            if shifted {
                                use popkern_core::geometry::curve::CurveEval;
                                let untorn = |orig: &NurbsCurve, work: &NurbsCurve, t: f64| {
                                    let o = orig.point_at(t);
                                    let w = work.point_at(t);
                                    let fold = (o.x + shift).rem_euclid(TWO_PI);
                                    let du = w.x - fold;
                                    (du - TWO_PI * (du / TWO_PI).round()).abs() < 1e-6
                                        && (w.y - o.y).abs() < 1e-6 * (1.0 + o.y.abs())
                                };
                                if !untorn(curve_uv, &working_curve, *t2) {
                                    if std::env::var_os("PK_XING_TRACE").is_some() {
                                        eprintln!(
                                            "      [xing-torn] he {he_handle:?} t2={t2:.4} rejected"
                                        );
                                    }
                                    continue;
                                }
                            }
                            if std::env::var_os("PK_XING_TRACE").is_some() {
                                use popkern_core::geometry::curve::CurveEval;
                                let q = working_curve.point_at(*t2);
                                let qb = boundary_pcurve.point_at(*t1);
                                eprintln!(
                                    "      [xing] he {he_handle:?} t2={t2:.4} uv=({:.4},{:.4}) bnd-uv=({:.4},{:.4})",
                                    q.x, q.y, qb.x, qb.y
                                );
                            }
                            found.push(UvCrossing {
                                t_curve: *t2,
                                boundary_he: *he_handle,
                                t_boundary: *t1,
                            });
                        }
                        CurveIntersection::Overlap { t1_range, t2_range } => {
                            // A COLLINEAR overlap with a boundary
                            // edge means the curve runs ALONG the
                            // face's rim for that stretch. The
                            // honest entry/exit candidates are the
                            // overlap ENDS — the old single midpoint
                            // "crossing" clipped a flush gusset
                            // cap's plate record to a fake corner
                            // 1.2 mm inside the true one and the
                            // union ring never closed.
                            found.push(UvCrossing {
                                t_curve: t2_range.0,
                                boundary_he: *he_handle,
                                t_boundary: t1_range.0,
                            });
                            found.push(UvCrossing {
                                t_curve: t2_range.1,
                                boundary_he: *he_handle,
                                t_boundary: t1_range.1,
                            });
                        }
                    }
                }
            }
        }
        found
    };

    let period_of = || -> f64 {
        use popkern_core::geometry::SurfaceEval as _;
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .and_then(|s| s.geometry.parameter_periodicity().0)
            .unwrap_or(std::f64::consts::TAU)
    };

    crossings.extend(collect_crossings(needs_shift));
    if crossings.is_empty() && surface_is_periodic {
        // The precise curve×curve intersector can miss real
        // crossings when the query's UV extent dwarfs the trim (a
        // near-axis-parallel plane×cylinder ellipse spans ±1000 in v
        // while the face is 40 wide — measured: the gusset-web weld
        // line found no crossings here, fell through to the seam
        // dance, and the torn-polygon crossings clipped the record
        // to the phantom re-crossing). Rescue with an in/out
        // TRANSITION SCAN against the face containment predicate —
        // honest regardless of pcurve period branches. Scoped to
        // queries whose UV extent dwarfs the trim (the starvation
        // signature): firing it on comparable-scale pairs changed
        // zero-crossing containment verdicts the cone suites are
        // tuned to.
        let bbox_diag = |pts: &mut dyn Iterator<Item = (f64, f64)>| -> f64 {
            let (mut lu, mut hu, mut lv, mut hv) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
            for (u, v) in pts {
                lu = lu.min(u);
                hu = hu.max(u);
                lv = lv.min(v);
                hv = hv.max(v);
            }
            if lu > hu {
                return 0.0;
            }
            ((hu - lu).powi(2) + (hv - lv).powi(2)).sqrt()
        };
        let q_diag = bbox_diag(&mut curve_uv.control_points.iter().map(|p| (p.x, p.y)));
        let b_diag = bbox_diag(
            &mut boundary_nurbs
                .iter()
                .flat_map(|(_, pc)| pc.control_points.iter().map(|p| (p.x, p.y))),
        );
        if q_diag.is_finite() && b_diag > 0.0 && q_diag > 8.0 * b_diag {
            crossings.extend(sampled_uv_crossings(
                brep,
                face,
                curve_uv,
                &boundary_nurbs,
                period_of(),
                tolerance,
            ));
        }
    }
    if crossings.is_empty() && surface_is_periodic {
        // The branch heuristic is a guess; when the chosen branch
        // finds NO crossings, a curve genuinely crossing the boundary
        // in the other branch would be misjudged by the containment
        // probe and clipped away whole (measured: two of the four
        // fillet-cylinder ellipses of a rotated lip union sat one
        // period from their quarter-arc charts and were discarded, so
        // the island ring never closed). Retry the other branch and
        // trust whichever finds real crossings — verified like the
        // first pass.
        crossings.extend(collect_crossings(!needs_shift));
    }

    crossings.sort_by(|a, b| {
        a.t_curve
            .partial_cmp(&b.t_curve)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Crossings must lie ON the curve: the pcurve-pair intersector
    // works on basis extensions, and hits beyond the curve's own
    // parameter window are not boundary crossings (measured: a fully
    // interior chord got one phantom crossing at t ≈ 4π and was
    // clipped away as Outside, losing socket segments on the second
    // boolean of a part).
    let (ct0, ct1) = curve_uv.parameter_range();
    let span = (ct1 - ct0).abs().max(1e-12);
    let tol_t = span * 1e-6;
    crossings.retain(|c| c.t_curve >= ct0.min(ct1) - tol_t && c.t_curve <= ct0.max(ct1) + tol_t);
    Some(crossings)
}

/// Extract a WRAPPED window `[t_enter, t_exit]` (with `t_exit`
/// beyond the curve's range end) from a CLOSED curve: the kept span
/// runs through the parameter seam. A closed curve's anchor is
/// arbitrary — extract `[t_enter, t_max]` and `[t_min, t_exit −
/// period]` and re-anchor the second piece after the first (C0
/// splice; the pieces meet AT the seam point). Wrapping is an
/// explicit caller decision (a seam-running [`UvClipResult::Windows`]
/// entry keeps its exit beyond the range end) — the plain extraction
/// below keeps its historical clamping for the out-of-range exits
/// ordinary callers pass.
pub(crate) fn clip_curve_wrapped(
    curve: &NurbsCurve,
    t_enter: f64,
    t_exit: f64,
) -> Option<NurbsCurve> {
    let (t_min, t_max) = curve.parameter_range();
    if t_exit <= t_max + 1e-12 {
        return clip_curve_to_crossings(curve, t_enter, t_exit);
    }
    if t_enter >= t_max - 1e-12 {
        return clip_curve_to_crossings(curve, t_min, t_exit - (t_max - t_min));
    }
    let period = t_max - t_min;
    let t_wrap = t_exit - period;
    if t_wrap <= t_min + 1e-12 {
        return clip_curve_to_crossings(curve, t_enter, t_max);
    }
    let head = clip_curve_to_crossings(curve, t_enter, t_max)?;
    let tail = clip_curve_to_crossings(curve, t_min, t_wrap)?;
    let join_tol = {
        use popkern_core::geometry::curve::CurveEval as _;
        let p0 = curve.point_at(t_min);
        let p1 = curve.point_at(t_max);
        p0.distance_to(&p1) + 1e-6
    };
    head.concatenate_clamped(&tail, join_tol).ok()
}

pub(crate) fn clip_curve_to_crossings(
    curve: &NurbsCurve,
    t_enter: f64,
    t_exit: f64,
) -> Option<NurbsCurve> {
    let (t_min, t_max) = curve.parameter_range();
    let t0 = t_enter.max(t_min);
    let t1 = t_exit.min(t_max);
    if t1 - t0 < 1e-15 {
        return None;
    }
    if curve.degree == 1 {
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);
        return Some(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p0, p1],
            weights: vec![1.0, 1.0],
        });
    }
    let eps = (t_max - t_min) * 1e-10;
    let right = if t0 - t_min > eps {
        curve.split_at(t0).ok()?.1
    } else {
        curve.clone()
    };
    let (_, r_max) = right.parameter_range();
    if t1 < r_max - eps {
        right.split_at(t1).ok().map(|(left, _)| left)
    } else {
        Some(right)
    }
}

pub(super) fn uv_range_from_crossings(
    crossings: &[UvCrossing],
    curve: &NurbsCurve,
    is_closed: bool,
) -> Option<(f64, f64)> {
    let (t_min, t_max) = curve.parameter_range();
    match crossings.len() {
        0 if is_closed => Some((t_min, t_max)),
        0 => None,
        1 => None,
        _ => {
            let t_enter = crossings[0].t_curve;
            let t_exit = crossings[crossings.len() - 1].t_curve;
            if t_exit - t_enter < 1e-15 {
                if is_closed {
                    return Some((t_min, t_max));
                }
                return None;
            }
            Some((t_enter, t_exit))
        }
    }
}

pub(super) enum UvClipResult {
    NoPcurves,
    Outside,
    /// In-face parameter windows, ascending and disjoint. On a CLOSED
    /// curve a window may run through the parameter seam: its exit
    /// then lies beyond the range end by design, and extraction
    /// re-anchors through the seam (`clip_curve_wrapped`).
    Windows(Vec<(f64, f64)>),
}

pub(super) fn uv_clip_one_face(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    is_closed: bool,
    open_windows_allowed: bool,
) -> UvClipResult {
    let trace = std::env::var_os("PK_CLIP_TRACE").is_some();
    match find_uv_crossings(brep, face, curve_uv, &brep.tolerance) {
        None => UvClipResult::NoPcurves,
        Some(crossings) if crossings.is_empty() => {
            if trace {
                use popkern_core::geometry::curve::CurveEval;
                let (t0, t1) = curve_uv.parameter_range();
                let mid = curve_uv.point_at(0.5 * (t0 + t1));
                let p0 = curve_uv.point_at(t0);
                let p1 = curve_uv.point_at(t1);
                eprintln!(
                    "    [uvclip-zero] {face:?} uv ({:.1},{:.1})->({:.1},{:.1}) mid=({:.1},{:.1}) inside={}",
                    p0.x,
                    p0.y,
                    p1.x,
                    p1.y,
                    mid.x,
                    mid.y,
                    crate::split_inner_loop::point_inside_face_uv(brep, face, mid, &brep.tolerance)
                );
                let ring = crate::clip::face_uv_boundary(brep, face);
                let mut lo = (f64::MAX, f64::MAX);
                let mut hi = (f64::MIN, f64::MIN);
                let mut n = 0;
                for (_, pc) in ring.iter().flatten() {
                    let (a, b) = pc.parameter_range();
                    for i in 0..4 {
                        let q = pc.point_at(a + (b - a) * i as f64 / 3.0);
                        lo = (lo.0.min(q.x), lo.1.min(q.y));
                        hi = (hi.0.max(q.x), hi.1.max(q.y));
                        n += 1;
                    }
                }
                eprintln!(
                    "    [uvclip-zero] {face:?} boundary uv-bbox ({:.1},{:.1})..({:.1},{:.1}) from {n} samples",
                    lo.0, lo.1, hi.0, hi.1
                );
            }
            // ZERO boundary crossings is ambiguous: the curve lies
            // entirely inside the face OR entirely outside it. The
            // historical rules — open ⇒ Outside, closed ⇒ inside —
            // both misfire (measured: fully-interior wall×cap
            // segments of a pocket cut into a cylinder cap were
            // discarded, while a plane×wall circle wholly outside
            // the tool's trimmed bottom face was kept whole).
            // Disambiguate by containment of the curve midpoint.
            use popkern_core::geometry::curve::CurveEval;
            let (t_min, t_max) = curve_uv.parameter_range();
            let mid = curve_uv.point_at(0.5 * (t_min + t_max));
            if crate::split_inner_loop::point_inside_face_uv(brep, face, mid, &brep.tolerance) {
                UvClipResult::Windows(vec![(t_min, t_max)])
            } else {
                UvClipResult::Outside
            }
        }
        Some(crossings) if crossings.len() == 1 && !is_closed => {
            // ONE boundary crossing: the curve END lies inside the
            // face — the input was pre-truncated by a padded
            // parameter box, so the in-face span is bounded by the
            // curve end, not by a second crossing (measured: a seat
            // tool wall's trim on a flipped-slope deck top was
            // truncated mid-face by the tool's own padded box; the
            // single surviving crossing hit the `1 => None` arm and
            // a REAL ring record was dropped as Outside). Decide
            // which side is inside by endpoint containment.
            use popkern_core::geometry::curve::CurveEval;
            let t_cross = crossings[0].t_curve;
            let (t_min, t_max) = curve_uv.parameter_range();
            let span = (t_max - t_min).max(1e-12);
            let inset = span * 1e-3;
            let start_in = crate::split_inner_loop::point_inside_face_uv(
                brep,
                face,
                curve_uv.point_at(t_min + inset),
                &brep.tolerance,
            );
            let end_in = crate::split_inner_loop::point_inside_face_uv(
                brep,
                face,
                curve_uv.point_at(t_max - inset),
                &brep.tolerance,
            );
            if trace {
                use popkern_core::geometry::curve::CurveEval as _;
                let p0 = curve_uv.point_at(t_min + inset);
                let p1 = curve_uv.point_at(t_max - inset);
                eprintln!(
                    "    [uvclip-1x] {face:?} t_cross={t_cross:.2} start=({:.2},{:.2})in={start_in} end=({:.2},{:.2})in={end_in}",
                    p0.x, p0.y, p1.x, p1.y
                );
            }
            match (start_in, end_in) {
                (true, false) => UvClipResult::Windows(vec![(t_min, t_cross)]),
                (false, true) => UvClipResult::Windows(vec![(t_cross, t_max)]),
                (false, false) => UvClipResult::Outside,
                // Both inside with one crossing = a numeric graze of
                // the boundary; keep the whole curve.
                (true, true) => UvClipResult::Windows(vec![(t_min, t_max)]),
            }
        }
        Some(crossings) => {
            // THREE OR MORE crossings on a CLOSED record: the curve
            // enters and leaves the face repeatedly — the crossing
            // hull below would span regions outside the face
            // (measured: a sphere×plane circle crossing a cube face's
            // bounds four times kept one arc running THROUGH the
            // out-of-face gap and over a foreign junction). Classify
            // every inter-crossing interval honestly; ambiguous
            // verdicts (boundary-hugging curves) fall back to the
            // tuned hull logic. OPEN records take windows only when
            // the caller MEASURED them safe (PK-114: transversality +
            // the boundary-hugging guard in decisive_windows) —
            // flush/coincident boundary-running traces keep the tuned
            // hull; splitting them starved the gap-substitute
            // completion (46 opens + 32% volume, flush gusset).
            if (is_closed || open_windows_allowed)
                && crossings.len() >= 3
                && let Some(ws) = super::clip_windows::decisive_windows(
                    brep, face, curve_uv, &crossings, is_closed,
                )
            {
                if trace {
                    eprintln!(
                        "    [uvclip-multi] {face:?} {} crossings -> {} decisive window(s)",
                        crossings.len(),
                        ws.len()
                    );
                }
                return if ws.is_empty() {
                    UvClipResult::Outside
                } else {
                    UvClipResult::Windows(ws)
                };
            }
            match uv_range_from_crossings(&crossings, curve_uv, is_closed) {
                Some((t_enter, t_exit)) => {
                    // CLOSED curves: the naive crossing hull can be the
                    // WRONG SIDE — when the inside span wraps through the
                    // parameter seam, the hull is exactly its complement
                    // (measured: a bore-mouth ellipse breaching a wall
                    // kept only the 8% sliver beyond the wall, t[0.46,
                    // 0.54], and discarded the 92% real mouth arc — the
                    // wall face's opening then assembled from the sliver
                    // and failed closure by 27.1). Probe the hull's
                    // midpoint: if it lies OUTSIDE the face, emit the
                    // complement as a WRAPPED window (exit beyond the
                    // range end; the extraction re-anchors through the
                    // seam).
                    use popkern_core::geometry::curve::CurveEval;
                    let (t_min, t_max) = curve_uv.parameter_range();
                    if is_closed && t_exit - t_enter > 1e-12 && (t_exit - t_enter) < (t_max - t_min)
                    {
                        // The complement is taken only on UNANIMOUS
                        // multi-sample evidence — every interior sample
                        // of the hull outside AND every interior sample
                        // of the complement inside. A single-midpoint
                        // test misfired on flush fixtures whose closed
                        // curves hug face boundaries: on-boundary
                        // containment verdicts are coin flips there, and
                        // one flipped midpoint inverted a correct hull
                        // (measured: 8 tuned union/chamfer gates broke).
                        // Boundary-hugging curves fail BOTH unanimity
                        // tests and keep the hull.
                        let period = t_max - t_min;
                        let inside = |t: f64| {
                            let mut tt = t;
                            if tt > t_max {
                                tt -= period;
                            }
                            crate::split_inner_loop::point_inside_face_uv(
                                brep,
                                face,
                                curve_uv.point_at(tt),
                                &brep.tolerance,
                            )
                        };
                        let samples = |a: f64, b: f64| -> [f64; 3] {
                            [a + (b - a) * 0.25, a + (b - a) * 0.5, a + (b - a) * 0.75]
                        };
                        let hull_all_out = samples(t_enter, t_exit).iter().all(|&t| !inside(t));
                        let comp_all_in =
                            samples(t_exit, t_enter + period).iter().all(|&t| inside(t));
                        if hull_all_out && comp_all_in {
                            if trace {
                                eprintln!(
                                    "    [uvclip-wrap] {face:?} hull ({t_enter:.3},{t_exit:.3}) outside — complement ({t_exit:.3},{:.3})",
                                    t_enter + period
                                );
                            }
                            return UvClipResult::Windows(vec![(t_exit, t_enter + period)]);
                        }
                    }
                    UvClipResult::Windows(vec![(t_enter, t_exit)])
                }
                None => UvClipResult::Outside,
            }
        }
    }
}
