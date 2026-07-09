use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
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
                            found.push(UvCrossing {
                                t_curve: *t2,
                                boundary_he: *he_handle,
                                t_boundary: *t1,
                            });
                        }
                        CurveIntersection::Overlap { t1_range, t2_range } => {
                            let t1_mid = (t1_range.0 + t1_range.1) * 0.5;
                            let t2_mid = (t2_range.0 + t2_range.1) * 0.5;
                            found.push(UvCrossing {
                                t_curve: t2_mid,
                                boundary_he: *he_handle,
                                t_boundary: t1_mid,
                            });
                        }
                    }
                }
            }
        }
        found
    };

    crossings.extend(collect_crossings(needs_shift));
    if crossings.is_empty() && surface_is_periodic {
        // The branch heuristic is a guess; when the chosen branch
        // finds NO crossings, a curve genuinely crossing the boundary
        // in the other branch would be misjudged by the containment
        // probe and clipped away whole (measured: two of the four
        // fillet-cylinder ellipses of a rotated lip union sat one
        // period from their quarter-arc charts and were discarded, so
        // the island ring never closed). Retry the other branch and
        // trust whichever finds real crossings.
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

enum UvClipResult {
    NoPcurves,
    Outside,
    Range(f64, f64),
}

fn uv_clip_one_face(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    is_closed: bool,
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
                UvClipResult::Range(t_min, t_max)
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
                (true, false) => UvClipResult::Range(t_min, t_cross),
                (false, true) => UvClipResult::Range(t_cross, t_max),
                (false, false) => UvClipResult::Outside,
                // Both inside with one crossing = a numeric graze of
                // the boundary; keep the whole curve.
                (true, true) => UvClipResult::Range(t_min, t_max),
            }
        }
        Some(crossings) => match uv_range_from_crossings(&crossings, curve_uv, is_closed) {
            Some((t_enter, t_exit)) => UvClipResult::Range(t_enter, t_exit),
            None => UvClipResult::Outside,
        },
    }
}

pub(crate) fn clip_curve_uv_space(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
) -> Option<CurveOnSurface> {
    if std::env::var_os("PK_CLIP_TRACE").is_some() {
        use popkern_core::geometry::curve::CurveEval;
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
    let clip_a = uv_clip_one_face(brep_a, face_a, &curve.curve_uv1, curve.is_closed);
    let clip_b = uv_clip_one_face(brep_b, face_b, &curve.curve_uv2, curve.is_closed);

    if matches!(clip_a, UvClipResult::Outside) || matches!(clip_b, UvClipResult::Outside) {
        if std::env::var_os("PK_CLIP_TRACE").is_some() {
            let name = |r: &UvClipResult| match r {
                UvClipResult::NoPcurves => "NoPcurves",
                UvClipResult::Outside => "Outside",
                UvClipResult::Range(_, _) => "Range",
            };
            eprintln!(
                "  [uvclip] dropped: A={} B={}",
                name(&clip_a),
                name(&clip_b)
            );
        }
        return None;
    }
    if matches!(clip_a, UvClipResult::NoPcurves) && matches!(clip_b, UvClipResult::NoPcurves) {
        return Some(curve.clone());
    }

    if std::env::var_os("PK_CLIP_TRACE").is_some() {
        let name = |r: &UvClipResult| match r {
            UvClipResult::NoPcurves => "NoPcurves".to_string(),
            UvClipResult::Outside => "Outside".to_string(),
            UvClipResult::Range(e, x) => format!("Range({e:.3},{x:.3})"),
        };
        eprintln!(
            "  [uvclip] {face_a:?}x{face_b:?}: A={} B={}",
            name(&clip_a),
            name(&clip_b)
        );
    }
    let (t_min, t_max) = curve.curve_3d.parameter_range();
    let (a_enter, a_exit) = match clip_a {
        UvClipResult::Range(enter, exit) => (enter, exit),
        _ => (t_min, t_max),
    };
    let (b_enter, b_exit) = match clip_b {
        UvClipResult::Range(enter, exit) => (enter, exit),
        _ => (t_min, t_max),
    };

    let t_enter = a_enter.max(b_enter);
    let t_exit = a_exit.min(b_exit);
    if t_exit - t_enter < 1e-15 {
        return None;
    }
    if curve.is_closed && (t_enter - t_min).abs() < 1e-12 && (t_exit - t_max).abs() < 1e-12 {
        return Some(curve.clone());
    }

    let clipped_3d = clip_curve_to_crossings(&curve.curve_3d, t_enter, t_exit)?;
    let clipped_uv1 = clip_curve_to_crossings(&curve.curve_uv1, t_enter, t_exit)?;
    let clipped_uv2 = clip_curve_to_crossings(&curve.curve_uv2, t_enter, t_exit)?;

    Some(CurveOnSurface {
        curve_3d: clipped_3d,
        curve_uv1: clipped_uv1,
        curve_uv2: clipped_uv2,
        is_closed: false,
    })
}
