use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

pub(super) const TWO_PI: f64 = std::f64::consts::TAU;

/// A crossing between a UV intersection curve and a face's UV boundary.
#[derive(Clone)]
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

/// Abstention band for boundary-riding containment samples: within
/// this uv distance of the true boundary contour, a sample is not
/// evidence for either side (record-fit slop puts genuinely
/// boundary-running curves this far off the contour). Chart-relative
/// — an absolute multiple of the model tolerance mis-scales across
/// charts (measured: 100×linear read 1.0 on a sphere chart in
/// radians and swallowed every sample of a mid-latitude ring).
fn uv_boundary_band(brep: &BRep, curve_uv: &NurbsCurve) -> f64 {
    let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for p in &curve_uv.control_points {
        lo_u = lo_u.min(p.x);
        hi_u = hi_u.max(p.x);
        lo_v = lo_v.min(p.y);
        hi_v = hi_v.max(p.y);
    }
    let extent = ((hi_u - lo_u).powi(2) + (hi_v - lo_v).powi(2)).sqrt();
    (extent * 1.0e-3).max(brep.tolerance.linear() * 10.0)
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
            // Disambiguate by containment — MAJORITY over five
            // interior samples, not the midpoint alone: a
            // near-tangent record skims a trim edge at one station
            // (the hex-corner chamfer cone's hyperbola apex grazes
            // the cone's r-flat rim exactly at the record midpoint),
            // and the single-sample verdict flipped Outside by the
            // boundary band's width while every other station sat
            // well inside — the sixth wall's corner cuts dropped
            // whole and the merge stranded 4 record sides (PK-464).
            // With zero crossings all samples lie on one side modulo
            // boundary noise, so the majority is the honest verdict.
            use popkern_core::geometry::curve::CurveEval;
            let (t_min, t_max) = curve_uv.parameter_range();
            // REGION containment (hole-aware): a record wholly inside
            // a hole is off the face (the vent-slot class); a record
            // crossing a hole has crossings and never reaches here.
            // BOUNDARY BAND + DENSE POLYGON (PK-465): the tuned
            // 8-per-pcurve polygon carries a chord sagitta far wider
            // than record-fit slop on a full-circle rim (0.2 uv on an
            // M6 tap's floor), so a genuine record hugging the rim
            // 0.05 inside voted Outside and the floor×base-wall
            // record dropped whole, stranding the hole-bottom lens.
            // Here the vote is rare (zero-crossing records only), so
            // it can afford honesty: samples within record-fit slop
            // of the TRUE boundary contour abstain (they are not
            // evidence for either side — a boundary-running record's
            // samples ALL abstain and the record drops, the PK-114
            // stance: the boundary already carries that topology);
            // the remaining samples vote against a boundary polygon
            // refined well under the band.
            let band = uv_boundary_band(brep, curve_uv);
            let boundary_ring = crate::clip::face_uv_boundary(brep, face);
            let boundary_d = |q: popkern_core::geometry::Point3| -> f64 {
                let mut best = f64::INFINITY;
                for (_, pc) in boundary_ring.iter().flatten() {
                    let (a, b) = pc.parameter_range();
                    if !(a.is_finite() && b.is_finite()) {
                        continue;
                    }
                    for i in 0..=64 {
                        let t = a + (b - a) * i as f64 / 64.0;
                        let s = pc.point_at(t);
                        let d = ((s.x - q.x).powi(2) + (s.y - q.y).powi(2)).sqrt();
                        if d < best {
                            best = d;
                        }
                    }
                }
                best
            };
            let mut inside_votes = 0usize;
            let mut counted = 0usize;
            for k in 1..=5 {
                let t = t_min + (t_max - t_min) * k as f64 / 6.0;
                let q = curve_uv.point_at(t);
                if boundary_d(q) < band {
                    continue;
                }
                counted += 1;
                if crate::split_inner_loop::point_inside_face_region_uv_with_density(
                    brep,
                    face,
                    q,
                    &brep.tolerance,
                    64,
                ) {
                    inside_votes += 1;
                }
            }
            if trace {
                eprintln!(
                    "    [uvclip-vote] {face:?} counted={counted} inside={inside_votes} band={band:.2e}"
                );
            }
            if counted == 0 {
                // Every sample rides the boundary band: the banded
                // vote knows nothing — fall back to the tuned
                // unbanded majority at the tuned density (declaring
                // Outside here dropped every record of a coarse-
                // tolerance sphere-cube subtraction, whose band floor
                // spans the whole cap, and the serve lied 6x on
                // volume — the tolerance-envelope pin).
                let mut inside = 0usize;
                for k in 1..=5 {
                    let t = t_min + (t_max - t_min) * k as f64 / 6.0;
                    if crate::split_inner_loop::point_inside_face_region_uv(
                        brep,
                        face,
                        curve_uv.point_at(t),
                        &brep.tolerance,
                    ) {
                        inside += 1;
                    }
                }
                if inside >= 3 {
                    UvClipResult::Windows(vec![(t_min, t_max)])
                } else {
                    UvClipResult::Outside
                }
            } else if 2 * inside_votes > counted {
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
            // DEGENERATE-HULL rescue (PK-465): an OPEN record whose
            // crossings all cluster at ONE boundary junction touched
            // the boundary tangentially — through the junction, not
            // across it — and the crossing hull collapses to a point,
            // discarding the real in-face span on one side (measured:
            // the angled tap's floor×flank record runs from a B-edge
            // junction 0.03 inside the floor rim, exits tangentially
            // THROUGH the rim junction, and its 2-crossing hull
            // (0.076,0.076) kept a zero-length stub — the hole-bottom
            // chain lost its closing link). The cluster discriminator
            // is GEOMETRIC — all crossings within the uv boundary
            // band of one point — because within record-fit slop of
            // the boundary two crossings are indistinguishable from
            // one junction by the same standard the containment
            // verdict below uses; a parametric epsilon (span×1e-6)
            // was a knife edge (measured: an unrelated registration
            // change moved this fixture's crossing pair from ~9e-7
            // to ~1.1e-6 apart in t and the rescue silently stopped
            // firing — the hole-bottom chain stranded again). A hull
            // narrower than the band is sub-slop: deciding it by
            // endpoint evidence is more honest than serving a
            // point-window stub. Sides are decided like the
            // 1-crossing arm, by endpoint containment — but with the
            // boundary-banded dense verdict, since these endpoints hug
            // the rim well inside the tuned polygon's sagitta. An
            // endpoint that abstains (rides the boundary) leaves the
            // tuned hull behavior untouched.
            {
                use popkern_core::geometry::curve::CurveEval;
                let (t_min, t_max) = curve_uv.parameter_range();
                let span = (t_max - t_min).max(1e-12);
                if !is_closed && crossings.len() >= 2 {
                    let t_first = crossings[0].t_curve;
                    let t_last = crossings[crossings.len() - 1].t_curve;
                    let band = uv_boundary_band(brep, curve_uv);
                    let q_first = curve_uv.point_at(t_first);
                    let clustered = crossings.iter().all(|c| {
                        let q = curve_uv.point_at(c.t_curve);
                        ((q.x - q_first.x).powi(2) + (q.y - q_first.y).powi(2)).sqrt() < band
                    });
                    if (t_last - t_first).abs() < span * 1e-6 || clustered {
                        let inset = span * 1e-3;
                        let boundary_ring = crate::clip::face_uv_boundary(brep, face);
                        let banded_inside = |q: popkern_core::geometry::Point3| -> Option<bool> {
                            let mut best = f64::INFINITY;
                            for (_, pc) in boundary_ring.iter().flatten() {
                                let (a, b) = pc.parameter_range();
                                if !(a.is_finite() && b.is_finite()) {
                                    continue;
                                }
                                for i in 0..=64 {
                                    let t = a + (b - a) * i as f64 / 64.0;
                                    let s = pc.point_at(t);
                                    let d = ((s.x - q.x).powi(2) + (s.y - q.y).powi(2)).sqrt();
                                    if d < best {
                                        best = d;
                                    }
                                }
                            }
                            if best < band {
                                return None;
                            }
                            Some(
                                crate::split_inner_loop::point_inside_face_region_uv_with_density(
                                    brep,
                                    face,
                                    q,
                                    &brep.tolerance,
                                    64,
                                ),
                            )
                        };
                        let start_in = banded_inside(curve_uv.point_at(t_min + inset));
                        let end_in = banded_inside(curve_uv.point_at(t_max - inset));
                        if trace {
                            eprintln!(
                                "    [uvclip-tangent] {face:?} t~{t_first:.3} start_in={start_in:?} end_in={end_in:?}"
                            );
                        }
                        match (start_in, end_in) {
                            (Some(true), Some(false)) => {
                                return UvClipResult::Windows(vec![(t_min, t_first)]);
                            }
                            (Some(false), Some(true)) => {
                                return UvClipResult::Windows(vec![(t_last, t_max)]);
                            }
                            (Some(true), Some(true)) => {
                                return UvClipResult::Windows(vec![(t_min, t_max)]);
                            }
                            (Some(false), Some(false)) => {
                                return UvClipResult::Outside;
                            }
                            // Abstaining endpoints: boundary-running —
                            // keep the tuned hull behavior below.
                            _ => {}
                        }
                    }
                }
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
