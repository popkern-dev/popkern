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
    let needs_shift = (is_near_seam(curve_uv)
        || boundary_nurbs.iter().any(|(_, pc)| is_near_seam(pc)))
        && (in_periodic_range(curve_uv)
            || boundary_nurbs.iter().any(|(_, pc)| in_periodic_range(pc)));
    let shift = std::f64::consts::PI;

    let working_curve = if needs_shift {
        shift_curve_u(curve_uv, shift)
    } else {
        curve_uv.clone()
    };

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

        let working_boundary = if needs_shift {
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
                        crossings.push(UvCrossing {
                            t_curve: *t2,
                            boundary_he: *he_handle,
                            t_boundary: *t1,
                        });
                    }
                    CurveIntersection::Overlap { t1_range, t2_range } => {
                        let t1_mid = (t1_range.0 + t1_range.1) * 0.5;
                        let t2_mid = (t2_range.0 + t2_range.1) * 0.5;
                        crossings.push(UvCrossing {
                            t_curve: t2_mid,
                            boundary_he: *he_handle,
                            t_boundary: t1_mid,
                        });
                    }
                }
            }
        }
    }

    crossings.sort_by(|a, b| {
        a.t_curve
            .partial_cmp(&b.t_curve)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

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
    match find_uv_crossings(brep, face, curve_uv, &brep.tolerance) {
        None => UvClipResult::NoPcurves,
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
    let clip_a = uv_clip_one_face(brep_a, face_a, &curve.curve_uv1, curve.is_closed);
    let clip_b = uv_clip_one_face(brep_b, face_b, &curve.curve_uv2, curve.is_closed);

    if matches!(clip_a, UvClipResult::Outside) || matches!(clip_b, UvClipResult::Outside) {
        return None;
    }
    if matches!(clip_a, UvClipResult::NoPcurves) && matches!(clip_b, UvClipResult::NoPcurves) {
        return Some(curve.clone());
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
