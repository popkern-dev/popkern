use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

pub(super) const TWO_PI: f64 = std::f64::consts::TAU;

/// A crossing between a UV intersection curve and a face's UV boundary.
pub(crate) struct UvCrossing {
    pub t_curve: f64,
    pub boundary_he: popkern_core::topology::HalfEdgeHandle,
    pub t_boundary: f64,
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

pub(crate) enum UvClipResult {
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
    let trace = popkern_core::trace::trace_enabled("clip");
    match super::clip_uv_crossings::find_uv_crossings(brep, face, curve_uv, &brep.tolerance) {
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
            // REGION containment (hole-aware): a record wholly inside
            // a hole is off the face (the vent-slot class); a record
            // crossing a hole has crossings and never reaches here.
            if crate::split_inner_loop::point_inside_face_region_uv(
                brep,
                face,
                mid,
                &brep.tolerance,
            ) {
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
                && let Some(ws) = super::clip_windows_decisive::decisive_windows(
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
                None => {
                    if trace {
                        eprintln!(
                            "    [uvclip-hull] {face:?} {} crossings, no hull -> Outside",
                            crossings.len()
                        );
                    }
                    UvClipResult::Outside
                }
            }
        }
    }
}
