//! The wire primitives (PK-1083): bodies of edges on no face, built in
//! the local frame and placed like the solids. A sweep or pipe takes a
//! wire's edge as its path, by name, without consuming the wire.
//!
//! Every curve is exact: lines and arcs analytic, the circle an arc of
//! a whole turn from the reference direction, the helix the sweep
//! kernel's exact helix, the spline the interpolant of its points at
//! chord-length parameters. Names (under the feature's operation): the
//! edges `Generated i` in order along the wire, the vertices
//! `Junction i` likewise.

use std::f64::consts::TAU;

use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::analytic_curves::circle_basis;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
use popkern_core::naming_criteria::{GeometricCriteria, capture_edge, capture_vertex};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;

use super::{failed, invalid};
use crate::error::ReplayError;

fn point(p: [f64; 3]) -> Point3 {
    Point3::new(p[0], p[1], p[2])
}

/// The arc about `center` in the plane of `normal` from `a` through a
/// turn of `sweep`.
fn arc(center: Point3, normal: Vec3, a: Point3, sweep: f64) -> CurveOnEdge {
    let (e1, e2) = circle_basis(normal);
    let r = a - center;
    let start = r.dot(&e2).patan2(r.dot(&e1));
    CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius: r.length(),
            start_angle: start,
            end_angle: start + sweep,
        },
        t_start: start,
        t_end: start + sweep,
    }
}

/// A chain of edges through `points` with `curves[i]` from point `i`
/// to point `i + 1` (wrapping when `closed`).
fn chain(points: &[Point3], curves: Vec<CurveOnEdge>, closed: bool, tolerance: Tolerance) -> BRep {
    let mut b = BRepBuilder::new(tolerance);
    let vs: Vec<_> = points.iter().map(|p| b.add_vertex(*p)).collect();
    for (i, curve) in curves.into_iter().enumerate() {
        let j = if closed { (i + 1) % vs.len() } else { i + 1 };
        let (eh, _, _) = b.add_edge(vs[i], vs[j]);
        b.set_edge_curve(eh, curve);
    }
    b.build_unchecked()
}

/// A polyline wire.
pub(super) fn polyline(
    shape: &str,
    points: &[[f64; 3]],
    closed: bool,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let min = if closed { 3 } else { 2 };
    if points.len() < min {
        return Err(invalid(
            shape,
            "points",
            points.len() as f64,
            &format!("needs at least {min} points"),
        ));
    }
    let pts: Vec<Point3> = points.iter().copied().map(point).collect();
    let n = if closed { pts.len() } else { pts.len() - 1 };
    let mut curves = Vec::with_capacity(n);
    for i in 0..n {
        let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
        let segment = CurveOnEdge::line_segment(a, b, tolerance.linear())
            .ok_or_else(|| invalid(shape, "points", i as f64, "two consecutive points coincide"))?;
        curves.push(segment);
    }
    Ok(chain(&pts, curves, closed, tolerance))
}

/// The arc through three points.
pub(super) fn three_point_arc(
    shape: &str,
    [start, via, end]: [[f64; 3]; 3],
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let (a, m, b) = (point(start), point(via), point(end));
    let (u, v) = (m - a, b - a);
    let normal = u.cross(&v);
    let n2 = normal.dot(&normal);
    if n2.sqrt() <= tolerance.linear() * (u.length() + v.length()) {
        return Err(invalid(
            shape,
            "via",
            0.0,
            "the three points are collinear: no arc passes through them",
        ));
    }
    // The circumcentre of the triangle a, m, b.
    let center = a + (normal.cross(&u) * v.dot(&v) + v.cross(&normal) * u.dot(&u)) * (0.5 / n2);
    let normal = normal * (1.0 / n2.sqrt());
    let (e1, e2) = circle_basis(normal);
    let angle = |p: Point3| {
        let r = p - center;
        r.dot(&e2).patan2(r.dot(&e1))
    };
    let sweep = (angle(b) - angle(a)).rem_euclid(TAU);
    Ok(chain(
        &[a, b],
        vec![arc(center, normal, a, sweep)],
        false,
        tolerance,
    ))
}

/// The whole circle about the local z from the reference direction.
pub(super) fn circle(radius: f64, tolerance: Tolerance) -> BRep {
    let start = Point3::new(radius, 0.0, 0.0);
    let curve = arc(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        start,
        TAU,
    );
    chain(&[start], vec![curve], true, tolerance)
}

/// A helix about the local z from the reference direction.
pub(super) fn helix(
    shape: &str,
    (radius, pitch, turns): (f64, f64, f64),
    left_handed: bool,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let spec = popkern_sweep::HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: radius,
        end_radius: radius,
        pitch: popkern_sweep::HelixPitch::Constant(pitch),
        turns,
        start_angle: 0.0,
        left_handed,
    };
    let (curve, _) =
        popkern_sweep::helix::build_helix_path(&spec).map_err(|e| failed(shape, "helix", e))?;
    Ok(nurbs_wire(curve, tolerance))
}

/// The spline interpolating `points` at chord-length parameters, or
/// the clamped uniform B-spline over them as control points.
pub(super) fn spline(
    shape: &str,
    points: &[[f64; 3]],
    degree: u32,
    control: bool,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let pts: Vec<Point3> = points.iter().copied().map(point).collect();
    if degree == 0 || pts.len() < degree as usize + 1 {
        return Err(invalid(
            shape,
            "points",
            pts.len() as f64,
            &format!(
                "a degree-{degree} spline needs at least {} points",
                degree + 1
            ),
        ));
    }
    if control {
        // Clamped uniform knots over the control polygon.
        let (n, p) = (pts.len(), degree as usize);
        let spans = n - p;
        let mut knots = vec![0.0; p + 1];
        knots.extend((1..spans).map(|i| i as f64 / spans as f64));
        knots.extend(vec![1.0; p + 1]);
        let curve = NurbsCurve {
            degree: p,
            knots,
            control_points: pts,
            weights: vec![1.0; n],
        };
        return Ok(nurbs_wire(curve, tolerance));
    }
    let mut params = vec![0.0];
    for w in pts.windows(2) {
        let step = w[0].distance_to(&w[1]);
        if step <= tolerance.linear() {
            return Err(invalid(
                shape,
                "points",
                0.0,
                "two consecutive points coincide",
            ));
        }
        params.push(params[params.len() - 1] + step);
    }
    let total = params[params.len() - 1];
    for t in &mut params {
        *t /= total;
    }
    let curve = popkern_core::geometry::fitting::interpolate_curve(&pts, &params, degree as usize)
        .map_err(|e| failed(shape, "spline", e))?;
    Ok(nurbs_wire(curve, tolerance))
}

/// One edge carrying `curve` between its end points.
fn nurbs_wire(curve: NurbsCurve, tolerance: Tolerance) -> BRep {
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = curve.parameter_range();
    let (a, b) = (curve.point_at(t0), curve.point_at(t1));
    let closed = a.distance_to(&b) <= tolerance.linear();
    let points: Vec<Point3> = if closed { vec![a] } else { vec![a, b] };
    let edge = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve),
        t_start: t0,
        t_end: t1,
    };
    chain(&points, vec![edge], closed, tolerance)
}

/// Name the wire's edges and vertices in order.
pub(super) fn name_wire(brep: &mut BRep, ctx: &mut TrackingContext) {
    let mut mint = |category: RoleCategory, index: usize| {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category,
                index: index as u32,
            },
            ancestors: Vec::new(),
        };
        ctx.record_generated(name.clone());
        name
    };
    let edges: Vec<_> = brep.edges.keys().collect();
    for (i, eh) in edges.into_iter().enumerate() {
        let name = mint(RoleCategory::Generated, i);
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let id = brep.name_pool.intern(name, criteria);
        brep.edges[eh].name = Some(id);
    }
    let vertices: Vec<_> = brep.vertices.keys().collect();
    for (i, vh) in vertices.into_iter().enumerate() {
        let name = mint(RoleCategory::Junction, i);
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let id = brep.name_pool.intern(name, criteria);
        brep.vertices[vh].name = Some(id);
    }
}
