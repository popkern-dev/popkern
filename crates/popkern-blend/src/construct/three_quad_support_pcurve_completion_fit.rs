//! Pcurve and edge-curve fitting helpers for cap/support seam completion.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::fitting::fit_curve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, nurbs_curve::NurbsCurve};

use super::support_pcurve_completion_quality::pcurve_surface_sample_gap;

pub(super) fn fit_nurbs_curve(points: &[Point3]) -> Option<NurbsCurve> {
    let params = (0..points.len())
        .map(|index| index as f64 / (points.len() - 1) as f64)
        .collect::<Vec<_>>();
    fit_curve(points, &params, 3, points.len()).ok()
}

pub(super) fn best_pcurve_for_surface_samples(
    surface: &SurfaceGeometry,
    uv_samples: &[Point3],
    target_samples: &[Point3],
) -> Pcurve {
    let linear = linear_pcurve_from_endpoints(uv_samples);
    if uv_samples_are_collinear(uv_samples) {
        return linear;
    }
    let linear_gap = pcurve_surface_sample_gap(surface, &linear, target_samples);
    let Some(fitted) = fit_nurbs_curve(uv_samples).map(Pcurve::from) else {
        return linear;
    };
    let fitted_gap = pcurve_surface_sample_gap(surface, &fitted, target_samples);
    if fitted_gap <= linear_gap {
        fitted
    } else {
        linear
    }
}

fn uv_samples_are_collinear(samples: &[Point3]) -> bool {
    if samples.len() < 3 {
        return true;
    }
    let start = samples[0];
    let end = samples[samples.len() - 1];
    let dx = end.x - start.x;
    let dy = end.y - start.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len <= 1.0e-12 {
        return samples.iter().all(|sample| {
            let sx = sample.x - start.x;
            let sy = sample.y - start.y;
            (sx * sx + sy * sy).sqrt() <= 1.0e-10
        });
    }
    samples.iter().all(|sample| {
        let sx = sample.x - start.x;
        let sy = sample.y - start.y;
        (sx * dy - sy * dx).abs() / len <= 1.0e-8
    })
}

fn linear_pcurve_from_endpoints(points: &[Point3]) -> Pcurve {
    let start = points[0];
    let end = points[points.len() - 1];
    Pcurve::linear((start.x, start.y), (end.x, end.y), 0.0, 1.0)
}

pub(super) fn line_geometry_from_endpoints(points: &[Point3]) -> CurveGeometry {
    let start = points[0];
    let end = points[points.len() - 1];
    CurveGeometry::Line {
        origin: start,
        direction: end - start,
    }
}
