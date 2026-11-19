//! The fillet overflow's sections and skin (PK-1189, PK-1192,
//! `strip_overflow`): the rolling strip's spine, the ball solved against
//! its two supports in each section plane, and exact sections (the ball's
//! rational arcs, or degree-5 G2 sections) skinned into one surface.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, PortableMath, Vec3};
use popkern_core::tolerance::Tolerance;

use crate::contact::{BallSide, ContactOptions, ContactRequest, solve_rolling_ball_contact};

/// The rolling strip's spine over a unit parameter: its tube-centre
/// circle (the arc its foot trim spans, or the whole ring when the trim
/// closes) or its axis over the trim's extent.
pub(super) enum FootPath {
    Ring {
        center: Point3,
        axis: Vec3,
        e1: Vec3,
        e2: Vec3,
        major: f64,
        theta0: f64,
        sweep: f64,
        closed: bool,
    },
    Axis {
        start: Point3,
        end: Point3,
    },
}

impl FootPath {
    pub(super) fn of(
        surface: &SurfaceGeometry,
        trim: &NurbsCurve,
        tolerance: &Tolerance,
    ) -> Result<Self, String> {
        let (t0, t1) = trim.parameter_range();
        let (a, z) = (trim.point_at(t0), trim.point_at(t1));
        match surface {
            SurfaceGeometry::Torus {
                center,
                axis,
                major_radius,
                ..
            } => FootPath::ring(*center, *axis, *major_radius, trim, tolerance),
            SurfaceGeometry::Cylinder { origin, axis, .. } => {
                let on_axis = |p: Point3| *origin + *axis * (p - *origin).dot(axis);
                Ok(FootPath::Axis {
                    start: on_axis(a),
                    end: on_axis(z),
                })
            }
            _ => Err("the rolling strip is not a torus or cylinder".into()),
        }
    }

    /// The spine circle about `axis` through `center` (in the spine's
    /// plane) of radius `major`, over the angles the foot trim spans.
    fn ring(
        center: Point3,
        axis: Vec3,
        major: f64,
        trim: &NurbsCurve,
        tolerance: &Tolerance,
    ) -> Result<Self, String> {
        let (t0, t1) = trim.parameter_range();
        let (a, z) = (trim.point_at(t0), trim.point_at(t1));
        let flat = |p: Point3| {
            let d = p - center;
            d - axis * d.dot(&axis)
        };
        let e1 = flat(a).normalized().ok_or("foot trim starts on the axis")?;
        let e2 = axis.cross(&e1);
        let closed = a.distance_to(&z) <= tolerance.linear().max(1e-9) * 10.0;
        let angle = |p: Point3| {
            let d = flat(p);
            d.dot(&e2).patan2(d.dot(&e1))
        };
        let sweep = if closed {
            std::f64::consts::TAU
        } else {
            // The trim's own way round, read at its middle.
            let mid = popkern_core::geometry::chart::ops::reduce(
                angle(trim.point_at(0.5 * (t0 + t1))),
                std::f64::consts::TAU,
            );
            let end =
                popkern_core::geometry::chart::ops::wrap_into(angle(z), 0.0, std::f64::consts::TAU);
            if mid <= end {
                end
            } else {
                end - std::f64::consts::TAU
            }
        };
        Ok(FootPath::Ring {
            center,
            axis,
            e1,
            e2,
            major,
            theta0: 0.0,
            sweep,
            closed,
        })
    }

    /// A swept strip's spine read off its edge: a circular edge's axis
    /// with the ball solved against `supports` (wall, face) in the
    /// meridian plane through the trim's start; a straight edge's line
    /// offset to the balls solved at the trim's ends.
    pub(super) fn of_edge(
        edge: &popkern_core::geometry::curve::CurveOnEdge,
        trim: &NurbsCurve,
        (wall, face): (&SurfaceGeometry, &SurfaceGeometry),
        r: f64,
        tolerance: &Tolerance,
    ) -> Result<Self, String> {
        use popkern_core::geometry::curve::CurveGeometry;
        let (t0, t1) = trim.parameter_range();
        let (a, z) = (trim.point_at(t0), trim.point_at(t1));
        let ball = |p: Point3, plane_normal: Vec3| -> Result<Point3, String> {
            let uv = (face as &dyn SurfaceEval)
                .closest_point(&p)
                .ok_or("no seed on F")?;
            let seed = p + face.normal_at(uv.0, uv.1) * r;
            let seed = if solve_section(wall, face, r, &seed, &plane_normal, tolerance).is_ok() {
                seed
            } else {
                p - face.normal_at(uv.0, uv.1) * r
            };
            Ok(solve_section(wall, face, r, &seed, &plane_normal, tolerance)?.center)
        };
        match &edge.geometry {
            CurveGeometry::Circle { center, normal, .. }
            | CurveGeometry::Arc { center, normal, .. } => {
                let axis = normal.normalized().ok_or("degenerate edge circle")?;
                let d = a - *center;
                let radial = (d - axis * d.dot(&axis))
                    .normalized()
                    .ok_or("trim on the axis")?;
                let c0 = ball(a, axis.cross(&radial))?;
                let h = (c0 - *center).dot(&axis);
                let dc = c0 - *center;
                let major = (dc - axis * h).length();
                FootPath::ring(*center + axis * h, axis, major, trim, tolerance)
            }
            CurveGeometry::Line { direction, .. } => {
                let dir = direction.normalized().ok_or("degenerate edge line")?;
                Ok(FootPath::Axis {
                    start: ball(a, dir)?,
                    end: ball(z, dir)?,
                })
            }
            _ => Err("the rolling strip's edge is neither a circle nor a line".into()),
        }
    }

    pub(super) fn spine(&self, u: f64) -> Point3 {
        match self {
            FootPath::Ring {
                center,
                e1,
                e2,
                major,
                theta0,
                sweep,
                ..
            } => {
                let th = theta0 + sweep * u;
                *center + (*e1 * th.pcos() + *e2 * th.psin()) * *major
            }
            FootPath::Axis { start, end } => *start + (*end - *start) * u,
        }
    }

    pub(super) fn closed(&self) -> bool {
        matches!(self, FootPath::Ring { closed: true, .. })
    }
}

/// The rolling strip's ball centre at its foot point `p` and its section
/// plane's normal: the spine circle's meridian plane, or the spine line's
/// cross-section plane.
pub(super) fn rolling_section(path: &FootPath, p: &Point3) -> Result<(Point3, Vec3), String> {
    match path {
        FootPath::Ring {
            center,
            axis,
            major,
            ..
        } => {
            let d = *p - *center;
            let radial = (d - *axis * d.dot(axis))
                .normalized()
                .ok_or("foot point on the strip's axis")?;
            let normal = axis
                .cross(&radial)
                .normalized()
                .ok_or("degenerate meridian")?;
            Ok((*center + radial * *major, normal))
        }
        FootPath::Axis { start, end } => {
            let dir = (*end - *start).normalized().ok_or("degenerate spine")?;
            Ok((*start + dir * (*p - *start).dot(&dir), dir))
        }
    }
}

/// A support's seed: its chart point nearest the spine, the ball's side
/// and the chart window searched.
pub(super) type Seed = ((f64, f64), BallSide, ((f64, f64), (f64, f64)));

/// The ball of radius `r` touching `wall` and `onto` in the plane through
/// `spine` with normal `plane_normal`, on the side of each surface its
/// seed centre `spine` lies.
pub(super) fn solve_section(
    wall: &SurfaceGeometry,
    onto: &SurfaceGeometry,
    r: f64,
    spine: &Point3,
    plane_normal: &Vec3,
    tolerance: &Tolerance,
) -> Result<crate::contact::ContactSolution, String> {
    let seed = |surface: &SurfaceGeometry| -> Result<Seed, String> {
        let uv = (surface as &dyn SurfaceEval)
            .closest_point(spine)
            .ok_or("no seed on a support")?;
        let foot = surface.point_at(uv.0, uv.1);
        let side = if (*spine - foot).dot(&surface.normal_at(uv.0, uv.1)) >= 0.0 {
            BallSide::Outside
        } else {
            BallSide::Inside
        };
        let (pu, pv) =
            popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(surface);
        let half_u = pu.map_or(4.0 * r, |p| 0.25 * p);
        let half_v = pv.map_or(4.0 * r, |p| 0.25 * p);
        Ok((
            uv,
            side,
            (
                (uv.0 - half_u, uv.0 + half_u),
                (uv.1 - half_v, uv.1 + half_v),
            ),
        ))
    };
    let (uv0, side_0, range0) = seed(wall)?;
    let (uv1, side_1, range1) = seed(onto)?;
    let request = ContactRequest {
        surface_0: wall,
        surface_1: onto,
        side_0,
        side_1,
        radius: r,
        initial_uv0: uv0,
        uv0_range: range0,
        initial_uv1: uv1,
        uv1_range: range1,
        spine_point: *spine,
        spine_plane_normal: *plane_normal,
    };
    let solved = solve_rolling_ball_contact(&request, &ContactOptions::default(), tolerance, None)
        .map_err(|e| format!("section contact did not converge: {e:?}"))?;
    if solved.residual > tolerance.linear().max(1e-9) {
        return Err(format!("section contact residual {:.3e}", solved.residual));
    }
    Ok(solved)
}

/// The ball's rational arc from its contact on the convex strip (v = 0)
/// to its contact on the wall (v = 1).
pub(super) fn arc_section(solved: &crate::contact::ContactSolution) -> Result<NurbsCurve, String> {
    let (cps, w) =
        crate::fillet_surface::rational_arc_row(&solved.center, &solved.contact1, &solved.contact0)
            .ok_or("no rational arc between a section's contacts")?;
    Ok(NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: cps.to_vec(),
        weights: vec![1.0, w, 1.0],
    })
}

/// The share of the convex strip's curvature in a G2 section at a
/// fraction `f` of the overflow: 0 at its ends (the face's, so the
/// section is the rolling strip's own), 1 from a fifth of the way in,
/// smooth (C1) between.
pub(super) fn g2_ramp(f: f64) -> f64 {
    const RAMP: f64 = 0.2;
    let x = (f.min(1.0 - f) / RAMP).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

/// The degree-5 G2 section of the ball `solved` (contact0 on the wall,
/// contact1 on the convex strip) from the strip (v = 0) to the wall
/// (v = 1), its curvature at the strip end the blend `share` of the
/// convex strip's and the face's (the two sections' control points
/// blended, interior ones linear in the end curvature).
pub(super) fn g2_section(
    solved: &crate::contact::ContactSolution,
    wall: &SurfaceGeometry,
    (face, onto): (&SurfaceGeometry, &SurfaceGeometry),
    share: f64,
) -> Result<NurbsCurve, String> {
    let against = |support: &SurfaceGeometry| -> Result<NurbsCurve, String> {
        let mut solution = solved.clone();
        solution.uv1 = (support as &dyn SurfaceEval)
            .closest_point(&solved.contact1)
            .ok_or("the contact does not project onto its support")?;
        let section =
            crate::g2_cross_section::g2_cross_section_bezier(&crate::g2_cross_section::G2Request {
                surface_0: wall,
                surface_1: support,
                solution: &solution,
            })
            .map_err(|e| format!("G2 section: {e:?}"))?;
        let mut curve = section.curve;
        curve.control_points.reverse();
        curve.weights.reverse();
        Ok(curve)
    };
    let to_onto = against(onto)?;
    if share >= 1.0 {
        return Ok(to_onto);
    }
    let to_face = against(face)?;
    let mut blended = to_onto.clone();
    for (b, (p, q)) in blended
        .control_points
        .iter_mut()
        .zip(to_face.control_points.iter().zip(&to_onto.control_points))
    {
        *b = *p + (*q - *p) * share;
    }
    Ok(blended)
}

/// Skin section curves sharing one degree and knot vector: each column of
/// homogeneous control points interpolated cubically across the sections
/// (one knot vector, so the columns stay one surface), the sections
/// themselves kept exact; `u` runs over the sections on [0, 1].
pub(super) fn skin_sections(sections: &[NurbsCurve]) -> Result<NurbsSurface, String> {
    let n = sections.len();
    let first = sections.first().ok_or("no sections")?;
    let num_v = first.control_points.len();
    if sections.iter().any(|c| {
        c.control_points.len() != num_v || c.knots != first.knots || c.degree != first.degree
    }) {
        return Err("sections differ in degree or knots".into());
    }
    let degree = 3.min(n - 1);
    let params: Vec<f64> = (0..n).map(|k| k as f64 / (n - 1) as f64).collect();
    let mut columns = Vec::with_capacity(num_v);
    for c in 0..num_v {
        let weighted: Vec<Point3> = sections
            .iter()
            .map(|s| {
                let (p, w) = (s.control_points[c], s.weights[c]);
                Point3::new(p.x * w, p.y * w, p.z * w)
            })
            .collect();
        let weights: Vec<Point3> = sections
            .iter()
            .map(|s| Point3::new(s.weights[c], 0.0, 0.0))
            .collect();
        let wc = popkern_core::geometry::fitting::interpolate_curve(&weighted, &params, degree)
            .map_err(|e| format!("skin column: {e}"))?;
        let ww = popkern_core::geometry::fitting::interpolate_curve(&weights, &params, degree)
            .map_err(|e| format!("skin weights: {e}"))?;
        columns.push((wc, ww));
    }
    let knots_u = columns[0].0.knots.clone();
    let num_u = columns[0].0.control_points.len();
    let mut control_points = Vec::with_capacity(num_u * num_v);
    let mut weights = Vec::with_capacity(num_u * num_v);
    for iu in 0..num_u {
        for (wc, ww) in &columns {
            let w = ww.control_points[iu].x;
            if w <= 0.0 {
                return Err("skin weight not positive".into());
            }
            let p = wc.control_points[iu];
            control_points.push(Point3::new(p.x / w, p.y / w, p.z / w));
            weights.push(w);
        }
    }
    Ok(NurbsSurface {
        degree_u: degree,
        degree_v: first.degree,
        knots_u,
        knots_v: first.knots.clone(),
        control_points,
        num_u,
        num_v,
        weights,
    })
}
