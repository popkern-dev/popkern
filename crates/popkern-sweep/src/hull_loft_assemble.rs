//! The hull loft's body, face by face (split from `hull_loft.rs`,
//! 500-line cap).

use super::*;

/// One section's edge between consecutive rulings.
#[derive(Clone, Copy)]
enum Run {
    /// The rulings share this section's end: no edge (a cone's apex, a
    /// triangle's tip).
    Point,
    /// A straight side.
    Line,
    /// An arc piece about `centre` of `radius` over the rulings' angles.
    Arc { centre: Point3, radius: f64 },
}

/// The edge of `section` between rulings at angles `from` and `to`
/// whose ends are `p` and `q`.
fn run_between(section: &Section, p: Point3, q: Point3, from: f64, to: f64, lin: f64) -> Run {
    if p.distance_to(&q) <= lin {
        return Run::Point;
    }
    if (to - from).abs() <= ANGLE_EPS {
        return Run::Line;
    }
    match section.arc_at(0.5 * (from + to)) {
        Some((centre, radius)) => Run::Arc { centre, radius },
        None => Run::Line,
    }
}

/// The run's curve from `p` to `q` (counter-clockwise about `normal`
/// for an arc).
fn run_curve(run: Run, p: Point3, q: Point3, normal: Vec3) -> Option<CurveOnEdge> {
    match run {
        Run::Point => None,
        Run::Line => Some(line(p, q)),
        Run::Arc { centre, radius } => Some(arc_between(
            &CircleFit {
                center: centre,
                normal,
                radius,
                max_deviation: 0.0,
                closed: false,
            },
            normal,
            p,
            q,
        )),
    }
}

/// The straight edge from `p` to `q`, parameter 0 to 1.
fn line(p: Point3, q: Point3) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p, q],
            weights: vec![1.0, 1.0],
        }),
        t_start: 0.0,
        t_end: 1.0,
    }
}

/// The rational arc about `centre` over angles `from`..`to`.
fn arc_nurbs(centre: Point3, radius: f64, normal: Vec3, from: f64, to: f64) -> NurbsCurve {
    CurveGeometry::Arc {
        center: centre,
        normal,
        radius,
        start_angle: from,
        end_angle: to,
    }
    .to_nurbs_in_range((from, to))
}

/// The ruled surface between the runs of a wall: rows paired control
/// point by control point (both arcs written over the same angles, or
/// an apex repeated with the arc's weights), `u` along the runs and `v`
/// from the first section to the second.
fn ruled_wall(
    a: (Run, Point3),
    b: (Run, Point3),
    normal: Vec3,
    (from, to): (f64, f64),
) -> Option<NurbsSurface> {
    let row = |run: Run, at: Point3| match run {
        Run::Arc { centre, radius } => Some(arc_nurbs(centre, radius, normal, from, to)),
        _ => {
            let _ = at;
            None
        }
    };
    let (ra, rb) = (row(a.0, a.1), row(b.0, b.1));
    let template = ra.clone().or_else(|| rb.clone())?;
    let point_row = |p: Point3| NurbsCurve {
        control_points: vec![p; template.control_points.len()],
        ..template.clone()
    };
    let ra = ra.unwrap_or_else(|| point_row(a.1));
    let rb = rb.unwrap_or_else(|| point_row(b.1));
    if ra.knots != rb.knots || ra.control_points.len() != rb.control_points.len() {
        return None;
    }
    let mut control_points = Vec::with_capacity(2 * ra.control_points.len());
    let mut weights = Vec::with_capacity(control_points.capacity());
    for i in 0..ra.control_points.len() {
        control_points.push(ra.control_points[i]);
        control_points.push(rb.control_points[i]);
        weights.push(ra.weights[i]);
        weights.push(rb.weights[i]);
    }
    Some(NurbsSurface {
        degree_u: ra.degree,
        degree_v: 1,
        knots_u: ra.knots.clone(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        num_u: ra.control_points.len(),
        num_v: 2,
        control_points,
        weights,
    })
}

/// The Newell normal of the closed polygon `pts`.
fn newell(pts: &[Point3]) -> Option<Vec3> {
    let n = pts.len();
    let mut acc = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let (p, q) = (pts[i], pts[(i + 1) % n]);
        acc = acc
            + Vec3::new(
                (p.y - q.y) * (p.z + q.z),
                (p.z - q.z) * (p.x + q.x),
                (p.x - q.x) * (p.y + q.y),
            );
    }
    acc.normalized()
}

/// The plane through `pts` facing `normal`.
fn plane(origin: Point3, normal: Vec3) -> SurfaceOnFace {
    let (u_axis, v_axis) = circle_basis(normal);
    SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        },
        orientation: true,
    }
}

/// The hull body through the cycle of `rulings`.
pub(super) fn assemble(
    rulings: &[(Point3, Point3, f64)],
    sa: &Section,
    sb: &Section,
    normal: Vec3,
    tolerance: Tolerance,
) -> Option<BRep> {
    let lin = tolerance.linear();
    let m = rulings.len();
    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);
    // Vertices, one per distinct point of each section along the cycle.
    let vertex = |pts: &mut Vec<(Point3, VertexHandle)>, p: Point3, b: &mut BRepBuilder| {
        if let Some(&(_, v)) = pts.iter().find(|(q, _)| q.distance_to(&p) <= lin) {
            return v;
        }
        let v = b.add_vertex(p);
        pts.push((p, v));
        v
    };
    let (mut va, mut vb) = (Vec::new(), Vec::new());
    let mut ends: Vec<(VertexHandle, VertexHandle)> = Vec::with_capacity(m);
    for &(a, b, _) in rulings {
        let x = vertex(&mut va, a, &mut builder);
        let y = vertex(&mut vb, b, &mut builder);
        ends.push((x, y));
    }
    // Rulings.
    let mut ruling: Vec<(HalfEdgeHandle, HalfEdgeHandle)> = Vec::with_capacity(m);
    for (k, &(a, b, _)) in rulings.iter().enumerate() {
        let (e, up, down) = builder.add_edge(ends[k].0, ends[k].1);
        builder.set_edge_curve(e, line(a, b));
        ruling.push((up, down));
    }
    // Section runs and walls.
    let span = |k: usize| {
        let (from, mut to) = (rulings[k].2, rulings[(k + 1) % m].2);
        if to < from - ANGLE_EPS {
            to += std::f64::consts::TAU;
        }
        (from, to)
    };
    let mut bottom: Vec<HalfEdgeHandle> = Vec::new();
    let mut top: Vec<HalfEdgeHandle> = Vec::new();
    for k in 0..m {
        let k1 = (k + 1) % m;
        let (from, to) = span(k);
        let (a0, b0, _) = rulings[k];
        let (a1, b1, _) = rulings[k1];
        let ra = run_between(sa, a0, a1, from, to, lin);
        let rb = run_between(sb, b0, b1, from, to, lin);
        let mut loop_ = vec![ruling[k].0];
        let mut corners = vec![a0, b0];
        if let Some(c) = run_curve(rb, b0, b1, normal) {
            let (e, fwd, back) = builder.add_edge(ends[k].1, ends[k1].1);
            corners.push(c.geometry.point_at(0.5 * (c.t_start + c.t_end)));
            builder.set_edge_curve(e, c);
            loop_.push(fwd);
            top.push(back);
        }
        loop_.push(ruling[k1].1);
        corners.extend([b1, a1]);
        if let Some(c) = run_curve(ra, a0, a1, normal) {
            let (e, fwd, back) = builder.add_edge(ends[k].0, ends[k1].0);
            corners.push(c.geometry.point_at(0.5 * (c.t_start + c.t_end)));
            builder.set_edge_curve(e, c);
            loop_.push(back);
            bottom.push(fwd);
        }
        if loop_.len() < 3 {
            return None;
        }
        let face = builder.add_face(shell);
        builder.add_loop(&loop_, face).ok()?;
        let facing = newell(&corners)?;
        let surface = match (ra, rb) {
            (Run::Arc { .. }, _) | (_, Run::Arc { .. }) => {
                let s = ruled_wall((ra, a0), (rb, b0), normal, (from, to))?;
                let ((u0, u1), _) = s.parameter_range();
                let n = s.normal_at(0.5 * (u0 + u1), 0.5).normalized()?;
                SurfaceOnFace {
                    geometry: SurfaceGeometry::Nurbs(s),
                    orientation: n.dot(&facing) > 0.0,
                }
            }
            _ => plane(a0, facing),
        };
        builder.set_face_surface(face, surface);
    }
    // Caps: the first section's runs in order, the second's reversed.
    top.reverse();
    for (runs, at, facing) in [(bottom, rulings[0].0, normal), (top, rulings[0].1, -normal)] {
        let face = builder.add_face(shell);
        builder.add_loop(&runs, face).ok()?;
        builder.set_face_surface(face, plane(at, facing));
    }
    let mut body = builder.build_unchecked();
    body.recompute_pcurves_where_missing().ok()?;
    crate::topology::outward::orient_outward(&mut body);
    body.validate_all().is_empty().then_some(body)
}
