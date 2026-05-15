//! Anchor orientation and subspan trimming for Gregory boundary edges.

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::gregory_boundary::reverse_nurbs;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, SurfaceSide};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;

pub(super) fn orient_boundary_to_anchors(
    boundary_curve: NurbsCurve,
    owner_pcurve: NurbsCurve,
    other_pcurve: NurbsCurve,
    anchors: &EdgeAnchors,
    owner_side: SurfaceSide,
    tol: &Tolerance,
) -> Option<Result<(NurbsCurve, NurbsCurve, NurbsCurve, f64), EdgeTripleError>> {
    let EdgeAnchors::Open { start, end } = anchors else {
        return None;
    };
    let (t_start, t_end) = boundary_curve.parameter_range();
    let tol_lin = tol.linear().max(1.0e-9);
    if let Some(oriented) = full_boundary_orientation(
        &boundary_curve,
        &owner_pcurve,
        &other_pcurve,
        *start,
        *end,
        tol_lin,
    ) {
        return Some(Ok(oriented));
    }
    Some(trim_boundary_to_anchor_subspan(TrimBoundaryRequest {
        boundary_curve,
        owner_pcurve,
        other_pcurve,
        start: *start,
        end: *end,
        owner_side,
        range: (t_start, t_end),
        tolerance: tol_lin,
    }))
}

fn full_boundary_orientation(
    boundary_curve: &NurbsCurve,
    owner_pcurve: &NurbsCurve,
    other_pcurve: &NurbsCurve,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> Option<(NurbsCurve, NurbsCurve, NurbsCurve, f64)> {
    let (t_start, t_end) = boundary_curve.parameter_range();
    let p_lo = boundary_curve.point_at(t_start);
    let p_hi = boundary_curve.point_at(t_end);
    let fwd_err = start.distance_to(&p_lo).max(end.distance_to(&p_hi));
    let rev_err = start.distance_to(&p_hi).max(end.distance_to(&p_lo));
    if fwd_err <= tolerance && fwd_err <= rev_err {
        Some((
            boundary_curve.clone(),
            owner_pcurve.clone(),
            other_pcurve.clone(),
            fwd_err,
        ))
    } else if rev_err <= tolerance {
        Some((
            reverse_nurbs(boundary_curve),
            reverse_nurbs(owner_pcurve),
            reverse_nurbs(other_pcurve),
            rev_err,
        ))
    } else {
        None
    }
}

struct TrimBoundaryRequest {
    boundary_curve: NurbsCurve,
    owner_pcurve: NurbsCurve,
    other_pcurve: NurbsCurve,
    start: Point3,
    end: Point3,
    owner_side: SurfaceSide,
    range: (f64, f64),
    tolerance: f64,
}

fn trim_boundary_to_anchor_subspan(
    request: TrimBoundaryRequest,
) -> Result<(NurbsCurve, NurbsCurve, NurbsCurve, f64), EdgeTripleError> {
    let TrimBoundaryRequest {
        boundary_curve,
        owner_pcurve,
        other_pcurve,
        start,
        end,
        owner_side,
        range,
        tolerance,
    } = request;
    let Some(start_param) = closest_curve_parameter(&boundary_curve, start, range.0, range.1)
    else {
        return Err(anchor_off_surface(owner_side, f64::INFINITY));
    };
    let Some(end_param) = closest_curve_parameter(&boundary_curve, end, range.0, range.1) else {
        return Err(anchor_off_surface(owner_side, f64::INFINITY));
    };
    let endpoint_err = boundary_curve
        .point_at(start_param)
        .distance_to(&start)
        .max(boundary_curve.point_at(end_param).distance_to(&end));
    if endpoint_err > tolerance || (end_param - start_param).abs() <= 1.0e-12 {
        return Err(anchor_off_surface(owner_side, endpoint_err));
    }
    let (lo, hi, reversed) = if start_param <= end_param {
        (start_param, end_param, false)
    } else {
        (end_param, start_param, true)
    };
    let curve = trim_nurbs_between(&boundary_curve, lo, hi)
        .ok_or_else(|| anchor_off_surface(owner_side, endpoint_err))?;
    let owner = trim_nurbs_between(&owner_pcurve, lo, hi)
        .ok_or_else(|| anchor_off_surface(owner_side, endpoint_err))?;
    let other = trim_nurbs_between(&other_pcurve, lo, hi)
        .ok_or_else(|| anchor_off_surface(owner_side, endpoint_err))?;
    if reversed {
        Ok((
            reverse_nurbs(&curve),
            reverse_nurbs(&owner),
            reverse_nurbs(&other),
            endpoint_err,
        ))
    } else {
        Ok((curve, owner, other, endpoint_err))
    }
}

fn anchor_off_surface(which: SurfaceSide, distance: f64) -> EdgeTripleError {
    EdgeTripleError::AnchorOffSurface { which, distance }
}

fn closest_curve_parameter(
    curve: &NurbsCurve,
    target: Point3,
    t_start: f64,
    t_end: f64,
) -> Option<f64> {
    let samples = 33usize;
    let mut best_index = 0usize;
    let mut best_distance = f64::INFINITY;
    for index in 0..samples {
        let f = index as f64 / (samples - 1) as f64;
        let t = t_start + (t_end - t_start) * f;
        let distance = curve.point_at(t).distance_to(&target);
        if distance < best_distance {
            best_distance = distance;
            best_index = index;
        }
    }
    if !best_distance.is_finite() {
        return None;
    }
    let step = (t_end - t_start) / (samples - 1) as f64;
    let mut lo = t_start + step * best_index.saturating_sub(1) as f64;
    let mut hi = t_start + step * (best_index + 1).min(samples - 1) as f64;
    for _ in 0..32 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if curve.point_at(m1).distance_to(&target) <= curve.point_at(m2).distance_to(&target) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    Some(0.5 * (lo + hi))
}

fn trim_nurbs_between(curve: &NurbsCurve, t_start: f64, t_end: f64) -> Option<NurbsCurve> {
    let (range_start, range_end) = curve.parameter_range();
    let eps = (range_end - range_start).abs().max(1.0) * 1.0e-12;
    let lo = t_start.clamp(range_start, range_end);
    let hi = t_end.clamp(range_start, range_end);
    if hi - lo <= eps {
        return None;
    }
    let mut trimmed = curve.clone();
    if hi < range_end - eps {
        trimmed = trimmed.split_at(hi).ok()?.0;
    }
    if lo > range_start + eps {
        trimmed = trimmed.split_at(lo).ok()?.1;
    }
    Some(trimmed)
}
