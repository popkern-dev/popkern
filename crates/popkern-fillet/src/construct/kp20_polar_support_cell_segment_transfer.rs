//! Exact Bezier segment transfer for KP20 support-cell boundary ownership.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;

use super::BoundarySide;

pub(super) fn exact_segment_transfer_candidate(
    patch_a: &RectangularGregoryBiNonic,
    side_a: BoundarySide,
    start_a: f64,
    end_a: f64,
    patch_b: &RectangularGregoryBiNonic,
    side_b: BoundarySide,
    start_b: f64,
    end_b: f64,
) -> Option<(RectangularGregoryBiNonic, RectangularGregoryBiNonic)> {
    let a_full = full_boundary_span(start_a, end_a);
    let b_full = full_boundary_span(start_b, end_b);
    if a_full == b_full {
        return None;
    }
    let mut candidate_a = patch_a.clone();
    let mut candidate_b = patch_b.clone();
    if a_full {
        let controls = trim_mate_subcurve_controls(patch_b, side_b, start_b, end_b, start_a, end_a);
        set_boundary_controls(&mut candidate_a, side_a, controls);
    } else {
        let controls = trim_mate_subcurve_controls(patch_a, side_a, start_a, end_a, start_b, end_b);
        set_boundary_controls(&mut candidate_b, side_b, controls);
    }
    Some((candidate_a, candidate_b))
}

fn full_boundary_span(a: f64, b: f64) -> bool {
    (close(a, 0.0) && close(b, 1.0)) || (close(a, 1.0) && close(b, 0.0))
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}

fn trim_mate_subcurve_controls(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
    mate_start: f64,
    mate_end: f64,
    full_start: f64,
    full_end: f64,
) -> [Point3; 10] {
    let row = boundary_controls(patch, side);
    let start = mate_parameter_for_full_parameter(mate_start, mate_end, full_start, full_end, 0.0);
    let end = mate_parameter_for_full_parameter(mate_start, mate_end, full_start, full_end, 1.0);
    bezier_subcurve(row, start, end)
}

fn mate_parameter_for_full_parameter(
    mate_start: f64,
    mate_end: f64,
    full_start: f64,
    full_end: f64,
    full_parameter: f64,
) -> f64 {
    let fraction = if (full_end - full_start).abs() <= 1.0e-12 {
        0.0
    } else {
        (full_parameter - full_start) / (full_end - full_start)
    };
    lerp_scalar(mate_start, mate_end, 1.0 - fraction)
}

fn boundary_controls(patch: &RectangularGregoryBiNonic, side: BoundarySide) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = match side {
            BoundarySide::UMin => patch.grid[0][index],
            BoundarySide::UMax => patch.grid[9][index],
            BoundarySide::VMin => patch.grid[index][0],
            BoundarySide::VMax => patch.grid[index][9],
        };
    }
    out
}

fn set_boundary_controls(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    controls: [Point3; 10],
) {
    for (index, point) in controls.into_iter().enumerate() {
        match side {
            BoundarySide::UMin => patch.grid[0][index] = point,
            BoundarySide::UMax => patch.grid[9][index] = point,
            BoundarySide::VMin => patch.grid[index][0] = point,
            BoundarySide::VMax => patch.grid[index][9] = point,
        }
    }
}

fn bezier_subcurve(controls: [Point3; 10], start: f64, end: f64) -> [Point3; 10] {
    if close(start, 0.0) && close(end, 1.0) {
        return controls;
    }
    if close(start, 1.0) && close(end, 0.0) {
        return reversed_controls(controls);
    }
    if end < start {
        return reversed_controls(bezier_subcurve(controls, end, start));
    }
    let (_, right) = split_bezier(controls, start.clamp(0.0, 1.0));
    let local_end = if close(start, 1.0) {
        0.0
    } else {
        ((end - start) / (1.0 - start)).clamp(0.0, 1.0)
    };
    let (segment, _) = split_bezier(right, local_end);
    segment
}

fn split_bezier(controls: [Point3; 10], t: f64) -> ([Point3; 10], [Point3; 10]) {
    let mut net = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    net[0] = controls;
    for level in 1..10 {
        for index in 0..10 - level {
            net[level][index] = lerp_point(net[level - 1][index], net[level - 1][index + 1], t);
        }
    }
    let mut left = [Point3::new(0.0, 0.0, 0.0); 10];
    let mut right = [Point3::new(0.0, 0.0, 0.0); 10];
    for index in 0..10 {
        left[index] = net[index][0];
        right[index] = net[9 - index][index];
    }
    (left, right)
}

fn reversed_controls(mut controls: [Point3; 10]) -> [Point3; 10] {
    controls.reverse();
    controls
}

fn lerp_point(a: Point3, b: Point3, t: f64) -> Point3 {
    a + (b - a) * t
}

fn lerp_scalar(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}
