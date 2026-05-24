//! NURBS export helpers for analytic pcurve variants.

use std::f64::consts::{PI, TAU};

use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::Pcurve;

pub(super) fn sampled_to_nurbs_approx(pcurve: &Pcurve) -> NurbsCurve {
    sampled_to_nurbs_with_points(pcurve, sampled_points(pcurve, 17))
}

pub(super) fn cylinder_arc_via_rational_quad_to_nurbs(pcurve: &Pcurve) -> NurbsCurve {
    let start = pcurve.point_at(pcurve.t_start);
    let mut end = pcurve.point_at(pcurve.t_end);
    unwrap_u_relative_to(start.x, &mut end.x);
    NurbsCurve {
        degree: 1,
        knots: vec![pcurve.t_start, pcurve.t_start, pcurve.t_end, pcurve.t_end],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    }
}

fn sampled_to_nurbs_with_points(pcurve: &Pcurve, control_points: Vec<Point3>) -> NurbsCurve {
    let params = (0..control_points.len())
        .map(|i| {
            let s = i as f64 / (control_points.len() as f64 - 1.0);
            pcurve.t_start + s * (pcurve.t_end - pcurve.t_start)
        })
        .collect::<Vec<_>>();
    match crate::geometry::fitting::interpolate_curve(&control_points, &params, 3) {
        Ok(c) => c,
        Err(_) => NurbsCurve {
            degree: 1,
            knots: vec![pcurve.t_start, pcurve.t_start, pcurve.t_end, pcurve.t_end],
            control_points: vec![
                control_points[0],
                control_points.last().copied().unwrap_or(control_points[0]),
            ],
            weights: vec![1.0, 1.0],
        },
    }
}

fn sampled_points(pcurve: &Pcurve, samples: usize) -> Vec<Point3> {
    let samples = samples.max(2);
    (0..samples)
        .map(|i| {
            let s = i as f64 / (samples as f64 - 1.0);
            let t = pcurve.t_start + s * (pcurve.t_end - pcurve.t_start);
            pcurve.point_at(t)
        })
        .collect()
}

fn unwrap_u_relative_to(reference_u: f64, u: &mut f64) {
    while *u - reference_u > PI {
        *u -= TAU;
    }
    while *u - reference_u < -PI {
        *u += TAU;
    }
}
