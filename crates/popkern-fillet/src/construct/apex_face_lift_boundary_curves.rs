//! Boundary-curve extraction for exact-boundary support-face lifts.

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

pub(super) fn compatible_unit_boundary_curves_for_support(
    brep: &BRep,
    support_face: FaceHandle,
    edge_a: EdgeHandle,
    edge_b: EdgeHandle,
) -> Result<(NurbsCurve, NurbsCurve), String> {
    let curve_a = support_lift_boundary_curve_as_unit_nurbs(brep, support_face, edge_a)?;
    let curve_b = support_lift_boundary_curve_as_unit_nurbs(brep, support_face, edge_b)?;
    let target_degree = curve_a.degree.max(curve_b.degree);
    let mut curve_a = elevate_boundary_curve(curve_a, target_degree)?;
    let mut curve_b = elevate_boundary_curve(curve_b, target_degree)?;
    let target_knots = merged_boundary_knots(&curve_a.knots, &curve_b.knots);
    insert_missing_boundary_knots(&mut curve_a, &target_knots)?;
    insert_missing_boundary_knots(&mut curve_b, &target_knots)?;
    Ok((curve_a, curve_b))
}

pub(super) fn write_edge_curve_from_boundary_row(
    brep: &mut BRep,
    edge: EdgeHandle,
    curve: NurbsCurve,
) {
    use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};

    let Some(edge_data) = brep.edges.get(edge) else {
        return;
    };
    let primary = edge_data.half_edge;
    let twin = brep.half_edges[primary].twin;
    let start = brep.vertices[brep.half_edges[primary].origin].position;
    let end = brep.vertices[brep.half_edges[twin].origin].position;
    let (range_start, range_end) = curve.parameter_range();
    let p_start = curve.point_at(range_start);
    let p_end = curve.point_at(range_end);
    let forward = start.distance_to(&p_start) + end.distance_to(&p_end)
        <= start.distance_to(&p_end) + end.distance_to(&p_start);
    let (t_start, t_end) = if forward {
        (range_start, range_end)
    } else {
        (range_end, range_start)
    };
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve),
        t_start,
        t_end,
    });
}

fn support_lift_boundary_curve_as_unit_nurbs(
    brep: &BRep,
    support_face: FaceHandle,
    edge: EdgeHandle,
) -> Result<NurbsCurve, String> {
    if let Some(curve) = cap_side_nurbs_boundary_curve(brep, support_face, edge) {
        return Ok(curve);
    }
    edge_curve_as_unit_nurbs(brep, edge)
}

fn cap_side_nurbs_boundary_curve(
    brep: &BRep,
    support_face: FaceHandle,
    edge: EdgeHandle,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::pcurve::PcurveGeometry;

    let edge_data = brep.edges.get(edge)?;
    let primary = edge_data.half_edge;
    let twin = brep.half_edges.get(primary)?.twin;
    let primary_face = brep.loops.get(brep.half_edges.get(primary)?.loop_)?.face;
    let half_edge = if primary_face == support_face {
        twin
    } else {
        primary
    };
    let cap_face = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face;
    if cap_face == support_face {
        return None;
    }
    let pcurve = brep.half_edges.get(half_edge)?.pcurve.as_ref()?;
    let surface = &brep.faces.get(cap_face)?.surface.as_ref()?.geometry;
    if let SurfaceGeometry::Nurbs(nurbs_surface) = surface
        && let PcurveGeometry::Linear { start, end } = &pcurve.geometry
        && let Some(curve) = nurbs_boundary_row_from_linear_pcurve(nurbs_surface, *start, *end)
    {
        return Some(curve);
    }
    sample_cap_side_boundary_curve(surface, pcurve)
}

fn sample_cap_side_boundary_curve(
    surface: &SurfaceGeometry,
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::fitting::interpolate_curve;

    let (t_start, t_end) = pcurve.parameter_range();
    if (t_end - t_start).abs() < 1.0e-12 {
        return None;
    }
    let samples = 25usize;
    let mut points = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let t = t_start + (t_end - t_start) * fraction;
        let uv = pcurve.point_at(t);
        points.push(surface.point_at(uv.x, uv.y));
        params.push(fraction);
    }
    interpolate_curve(&points, &params, 7).ok()
}

fn nurbs_boundary_row_from_linear_pcurve(
    surface: &NurbsSurface,
    start: (f64, f64),
    end: (f64, f64),
) -> Option<NurbsCurve> {
    let same = |a: f64, b: f64| (a - b).abs() <= 1.0e-8;
    let reversed = (end.0 + end.1) < (start.0 + start.1);
    let mut curve = if same(start.1, end.1) && same(start.1, 0.0) {
        nurbs_surface_u_boundary(surface, 0)
    } else if same(start.1, end.1) && same(start.1, 1.0) {
        nurbs_surface_u_boundary(surface, surface.num_v.checked_sub(1)?)
    } else if same(start.0, end.0) && same(start.0, 0.0) {
        nurbs_surface_v_boundary(surface, 0)
    } else if same(start.0, end.0) && same(start.0, 1.0) {
        nurbs_surface_v_boundary(surface, surface.num_u.checked_sub(1)?)
    } else {
        return None;
    };
    if reversed {
        reverse_unit_nurbs_curve(&mut curve);
    }
    Some(curve)
}

fn nurbs_surface_u_boundary(surface: &NurbsSurface, i_v: usize) -> NurbsCurve {
    let mut control_points = Vec::with_capacity(surface.num_u);
    let mut weights = Vec::with_capacity(surface.num_u);
    for i_u in 0..surface.num_u {
        let index = i_u * surface.num_v + i_v;
        control_points.push(surface.control_points[index]);
        weights.push(surface.weights[index]);
    }
    NurbsCurve {
        degree: surface.degree_u,
        knots: surface.knots_u.clone(),
        control_points,
        weights,
    }
}

fn nurbs_surface_v_boundary(surface: &NurbsSurface, i_u: usize) -> NurbsCurve {
    let mut control_points = Vec::with_capacity(surface.num_v);
    let mut weights = Vec::with_capacity(surface.num_v);
    for i_v in 0..surface.num_v {
        let index = i_u * surface.num_v + i_v;
        control_points.push(surface.control_points[index]);
        weights.push(surface.weights[index]);
    }
    NurbsCurve {
        degree: surface.degree_v,
        knots: surface.knots_v.clone(),
        control_points,
        weights,
    }
}

fn reverse_unit_nurbs_curve(curve: &mut NurbsCurve) {
    let (range_start, range_end) = {
        use popkern_core::geometry::curve::CurveEval;
        curve.parameter_range()
    };
    let sum = range_start + range_end;
    curve.control_points.reverse();
    curve.weights.reverse();
    curve.knots = curve.knots.iter().rev().map(|knot| sum - *knot).collect();
}

fn edge_curve_as_unit_nurbs(brep: &BRep, edge: EdgeHandle) -> Result<NurbsCurve, String> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::curve::CurveGeometry;

    let coe = brep
        .edges
        .get(edge)
        .and_then(|edge_data| edge_data.curve.as_ref())
        .ok_or_else(|| "boundary edge has no curve3d".to_string())?;
    if (coe.t_end - coe.t_start).abs() < 1.0e-12 {
        return Err("boundary edge has zero parameter span".into());
    }
    if let CurveGeometry::Hyperbola {
        center,
        semi_major,
        semi_minor,
        ..
    } = &coe.geometry
    {
        return exact_hyperbola_segment_as_unit_nurbs(
            *center,
            *semi_major,
            *semi_minor,
            coe.t_start,
            coe.t_end,
        );
    }

    let mut nurbs = coe.geometry.to_nurbs_in_range((coe.t_start, coe.t_end));
    nurbs = trim_boundary_nurbs_to_edge_span(nurbs, coe.t_start, coe.t_end)?;
    let (k_min, k_max) = nurbs.parameter_range();
    let span = k_max - k_min;
    if span.abs() < 1.0e-12 {
        return Err("boundary NURBS has zero parameter span".into());
    }
    for knot in &mut nurbs.knots {
        *knot = (*knot - k_min) / span;
    }
    Ok(nurbs)
}

fn exact_hyperbola_segment_as_unit_nurbs(
    center: Point3,
    semi_major: popkern_core::geometry::Vec3,
    semi_minor: popkern_core::geometry::Vec3,
    t_start: f64,
    t_end: f64,
) -> Result<NurbsCurve, String> {
    let u0 = (0.5 * t_start).tanh();
    let u1 = (0.5 * t_end).tanh();
    let a = u0;
    let b = u1 - u0;
    if b.abs() < 1.0e-12 {
        return Err("hyperbola boundary span is degenerate".into());
    }

    let c = [center.x, center.y, center.z];
    let maj = [semi_major.x, semi_major.y, semi_major.z];
    let min = [semi_minor.x, semi_minor.y, semi_minor.z];

    let w0 = 1.0 - a * a;
    let w1 = -2.0 * a * b;
    let w2 = -(b * b);
    let q0 = [
        c[0] * (1.0 - a * a) + maj[0] * (1.0 + a * a) + min[0] * (2.0 * a),
        c[1] * (1.0 - a * a) + maj[1] * (1.0 + a * a) + min[1] * (2.0 * a),
        c[2] * (1.0 - a * a) + maj[2] * (1.0 + a * a) + min[2] * (2.0 * a),
    ];
    let q1 = [
        c[0] * (-2.0 * a * b) + maj[0] * (2.0 * a * b) + min[0] * (2.0 * b),
        c[1] * (-2.0 * a * b) + maj[1] * (2.0 * a * b) + min[1] * (2.0 * b),
        c[2] * (-2.0 * a * b) + maj[2] * (2.0 * a * b) + min[2] * (2.0 * b),
    ];
    let q2 = [
        c[0] * (-(b * b)) + maj[0] * (b * b),
        c[1] * (-(b * b)) + maj[1] * (b * b),
        c[2] * (-(b * b)) + maj[2] * (b * b),
    ];

    let homogeneous = [
        (q0, w0),
        (
            [
                q0[0] + 0.5 * q1[0],
                q0[1] + 0.5 * q1[1],
                q0[2] + 0.5 * q1[2],
            ],
            w0 + 0.5 * w1,
        ),
        (
            [
                q0[0] + q1[0] + q2[0],
                q0[1] + q1[1] + q2[1],
                q0[2] + q1[2] + q2[2],
            ],
            w0 + w1 + w2,
        ),
    ];

    let mut control_points = Vec::with_capacity(3);
    let mut weights = Vec::with_capacity(3);
    for (xyz_w, weight) in homogeneous {
        if weight.abs() < 1.0e-12 {
            return Err("hyperbola boundary has near-zero rational weight".into());
        }
        control_points.push(Point3::new(
            xyz_w[0] / weight,
            xyz_w[1] / weight,
            xyz_w[2] / weight,
        ));
        weights.push(weight);
    }

    Ok(NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points,
        weights,
    })
}

fn trim_boundary_nurbs_to_edge_span(
    mut curve: NurbsCurve,
    t_start: f64,
    t_end: f64,
) -> Result<NurbsCurve, String> {
    use popkern_core::geometry::curve::CurveEval;

    if t_end < t_start {
        return Err("reversed boundary NURBS spans are not supported".into());
    }
    let (range_start, range_end) = curve.parameter_range();
    let eps = (range_end - range_start).abs().max(1.0) * 1.0e-12;
    let lo = t_start.clamp(range_start, range_end);
    let hi = t_end.clamp(range_start, range_end);
    if hi - lo <= eps {
        return Err("boundary NURBS trim span is degenerate".into());
    }
    if hi < range_end - eps {
        curve = curve
            .split_at(hi)
            .map_err(|err| format!("boundary upper trim failed: {err}"))?
            .0;
    }
    if lo > range_start + eps {
        curve = curve
            .split_at(lo)
            .map_err(|err| format!("boundary lower trim failed: {err}"))?
            .1;
    }
    Ok(curve)
}

fn elevate_boundary_curve(curve: NurbsCurve, target_degree: usize) -> Result<NurbsCurve, String> {
    if curve.degree == target_degree {
        Ok(curve)
    } else {
        curve
            .elevate_degree_to(target_degree)
            .map_err(|err| format!("boundary degree elevation failed: {err}"))
    }
}

fn merged_boundary_knots(knots_a: &[f64], knots_b: &[f64]) -> Vec<f64> {
    let mut unique = Vec::new();
    for &knot in knots_a.iter().chain(knots_b.iter()) {
        if !unique
            .iter()
            .any(|existing: &f64| (*existing - knot).abs() <= 1.0e-9)
        {
            unique.push(knot);
        }
    }
    unique.sort_by(|a, b| a.total_cmp(b));

    let mut merged = Vec::new();
    for knot in unique {
        let count_a = knot_multiplicity(knots_a, knot);
        let count_b = knot_multiplicity(knots_b, knot);
        merged.extend(std::iter::repeat_n(knot, count_a.max(count_b)));
    }
    merged
}

fn insert_missing_boundary_knots(
    curve: &mut NurbsCurve,
    target_knots: &[f64],
) -> Result<(), String> {
    let mut unique = Vec::new();
    for &knot in target_knots {
        if knot <= 1.0e-9 || knot >= 1.0 - 1.0e-9 {
            continue;
        }
        if unique
            .iter()
            .any(|existing: &f64| (*existing - knot).abs() <= 1.0e-9)
        {
            continue;
        }
        unique.push(knot);
    }

    let mut to_insert = Vec::new();
    for knot in unique {
        let target_count = knot_multiplicity(target_knots, knot);
        let current_count = knot_multiplicity(&curve.knots, knot);
        if target_count > current_count {
            to_insert.extend(std::iter::repeat_n(knot, target_count - current_count));
        }
    }
    if to_insert.is_empty() {
        return Ok(());
    }
    to_insert.sort_by(|a, b| a.total_cmp(b));
    *curve = curve
        .insert_knots(&to_insert)
        .map_err(|err| format!("boundary knot insertion failed: {err}"))?;
    Ok(())
}

fn knot_multiplicity(knots: &[f64], knot: f64) -> usize {
    knots
        .iter()
        .filter(|&&candidate| (candidate - knot).abs() <= 1.0e-9)
        .count()
}
