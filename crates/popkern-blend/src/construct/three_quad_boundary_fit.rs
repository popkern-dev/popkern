//! Degree-9 boundary fitting for K-P rectangular cap promotion.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::fitting::fit_curve;

use crate::blend::boundary::{BlendBoundary, BoundarySegment};

pub(super) fn fit_hex_boundary_cps_deg9(
    boundary: &BlendBoundary,
) -> Option<Vec<([Point3; 10], [Point3; 10])>> {
    if !boundary.segments.len().is_multiple_of(2) {
        return None;
    }
    let n = boundary.segments.len() / 2;
    let mut fitted = Vec::with_capacity(n);
    for k in 0..n {
        let row0 = fit_segment_deg9(boundary.segments.get(2 * k)?)?;
        let col9 = fit_segment_deg9(boundary.segments.get(2 * k + 1)?)?;
        fitted.push((row0, col9));
    }
    Some(fitted)
}

fn fit_segment_deg9(segment: &BoundarySegment) -> Option<[Point3; 10]> {
    let curve = match segment {
        BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => curve,
    };
    let (t0, t1) = curve.parameter_range();
    if (t1 - t0).abs() < 1.0e-15 {
        return None;
    }
    let mut samples = Vec::with_capacity(10);
    let mut params = Vec::with_capacity(10);
    for i in 0..10 {
        let s = (i as f64) / 9.0;
        samples.push(curve.point_at(t0 + s * (t1 - t0)));
        params.push(s);
    }
    let fitted = fit_curve(&samples, &params, 9, 10).ok()?;
    if fitted.control_points.len() != 10 {
        return None;
    }
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (i, cp) in fitted.control_points.iter().enumerate() {
        out[i] = *cp;
    }
    Some(out)
}
