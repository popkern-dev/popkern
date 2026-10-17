//! Trig-pair sphere solvers on a sphere whose chart pole is not world +Z
//! (PK-858).
//!
//! The sphere solvers write pcurves in the world-+Z chart, so PK-857 made
//! a turned sphere decline to the NURBS-via-`intersect_surfaces` path
//! (fitted pcurves). Here the pair is solved in the sphere's own frame:
//! a rigid rotation `T` about the sphere's centre takes its chart frame
//! (`orthonormal_basis(axis)`, axis) onto the world chart frame
//! (`orthonormal_basis(Z)`, Z), so the turned sphere's chart IS the world
//! chart of its image — its pcurve carries over unchanged. The partner
//! surface rotates with it; its pcurve re-anchors by the constant U shift
//! `BRep::transform` applies to a rotated rotational chart (a plane's
//! chart rotates exactly). The 3-D curve maps back by `T⁻¹`; a circle or
//! arc re-parameterises with its derived basis, and both pcurves follow
//! that parameter shift so the parameter map stays the identity.
//!
//! Declines (`None`, the previous NURBS fallback) when two spheres have
//! different poles (no single frame makes both world) or when the frames
//! differ in handedness (the rotation would be a reflection).

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve::CurveOnEdge;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, EdgeHint};
use crate::geometry::edge_triple::types::EdgeTriple;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::transform::Transform3;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// The sphere's chart pole, if `surface` is a sphere.
fn sphere_pole(surface: &SurfaceGeometry) -> Option<(Point3, Vec3)> {
    match surface {
        SurfaceGeometry::Sphere { center, axis, .. } => Some((*center, *axis)),
        _ => None,
    }
}

/// The rigid rotation about `center` taking (e1, e2, axis) onto
/// (f1, f2, Z), or `None` for a reflection.
fn frame_rotation(center: Point3, axis: Vec3) -> Option<Transform3> {
    let a = axis.normalized()?;
    let (e1, e2) = orthonormal_basis(a);
    let z = Vec3::Z;
    let (f1, f2) = orthonormal_basis(z);
    if e1.cross(&e2).dot(&a) * f1.cross(&f2).dot(&z) <= 0.0 {
        return None;
    }
    let (from, to) = ([e1, e2, a], [f1, f2, z]);
    let mut r = [[0.0f64; 3]; 3];
    for (i, row) in r.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = (0..3)
                .map(|k| {
                    let comp = |v: Vec3, n: usize| [v.x, v.y, v.z][n];
                    comp(to[k], i) * comp(from[k], j)
                })
                .sum();
        }
    }
    let c = [center.x, center.y, center.z];
    let t: Vec<f64> = (0..3)
        .map(|i| c[i] - (0..3).map(|j| r[i][j] * c[j]).sum::<f64>())
        .collect();
    Some(Transform3 {
        matrix: [
            [r[0][0], r[0][1], r[0][2], t[0]],
            [r[1][0], r[1][1], r[1][2], t[1]],
            [r[2][0], r[2][1], r[2][2], t[2]],
            [0.0, 0.0, 0.0, 1.0],
        ],
    })
}

/// The U shift a rotational chart's pcurves take under `t` (as
/// `BRep::transform` re-anchors them); 0 for other charts.
fn chart_u_shift(surface: &SurfaceGeometry, t: &Transform3) -> f64 {
    let axis = match surface {
        SurfaceGeometry::Cylinder { axis, .. }
        | SurfaceGeometry::Cone { axis, .. }
        | SurfaceGeometry::Sphere { axis, .. }
        | SurfaceGeometry::Torus { axis, .. } => *axis,
        _ => return 0.0,
    };
    let (e1, _) = orthonormal_basis(axis);
    let (f1, f2) = orthonormal_basis(t.apply_direction(&axis));
    let image = t.apply_direction(&e1);
    image.dot(&f2).patan2(image.dot(&f1))
}

fn map_anchors(anchors: &EdgeAnchors, t: &Transform3) -> EdgeAnchors {
    match anchors {
        EdgeAnchors::Open { start, end } => EdgeAnchors::Open {
            start: t.apply_point(start),
            end: t.apply_point(end),
        },
        EdgeAnchors::Closed { sample, direction } => EdgeAnchors::Closed {
            sample: t.apply_point(sample),
            direction: t.apply_direction(direction),
        },
        EdgeAnchors::Seam { sample, seam_axis } => EdgeAnchors::Seam {
            sample: t.apply_point(sample),
            seam_axis: *seam_axis,
        },
    }
}

fn map_hint(hint: &EdgeHint, t: &Transform3) -> EdgeHint {
    match hint {
        EdgeHint::Curve3D(curve) => {
            let mut curve = curve.clone();
            for p in &mut curve.control_points {
                *p = t.apply_point(p);
            }
            EdgeHint::Curve3D(curve)
        }
        other => other.clone(),
    }
}

/// The trig dispatch for a pair holding a turned sphere, solved in the
/// sphere's frame (see the module docs); `None` when it does not apply.
pub(super) fn try_in_sphere_frame(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let (center, axis) = match (sphere_pole(surf_a), sphere_pole(surf_b)) {
        (Some((ca, aa)), Some((_, ab))) => {
            // One frame must make both poles world: equal poles only.
            if aa.cross(&ab).length() > 1e-12 || aa.dot(&ab) <= 0.0 {
                return None;
            }
            (ca, aa)
        }
        (Some(p), None) | (None, Some(p)) => p,
        (None, None) => return None,
    };
    let t = frame_rotation(center, axis)?;
    let back = t.inverse().ok()?;
    // `T` takes the pole onto Z exactly in exact arithmetic; its image
    // lands an ulp off, and the dispatch tests the world chart bit for
    // bit — an unsnapped pole re-entered this frame path until the
    // stack overflowed (the rounded tetrahedron's rolling-ball corners).
    let moved = |s: &SurfaceGeometry| match crate::topology::brep::transform_surface(s, &t, 1.0) {
        SurfaceGeometry::Sphere { center, radius, .. } => SurfaceGeometry::Sphere {
            center,
            radius,
            axis: Vec3::Z,
        },
        other => other,
    };
    let (ta, tb) = (moved(surf_a), moved(surf_b));
    if !(super::analytic_trig::sphere_chart_is_world(&ta)
        && super::analytic_trig::sphere_chart_is_world(&tb))
    {
        return None;
    }
    let original = anchors.clone();
    let anchors = map_anchors(anchors, &t);
    let hint = hint.map(|h| map_hint(h, &t));
    let solved = super::analytic_trig::try_dispatch_trig(&ta, &tb, &anchors, hint.as_ref(), tol)?;
    Some(solved.map(|mut triple| {
        let mut curve = CurveOnEdge {
            geometry: triple.curve_3d.clone(),
            t_start: triple.t_start,
            t_end: triple.t_end,
        };
        crate::topology::brep::transform_curve_on_edge(&mut curve, &back, 1.0);
        let shift = curve.t_start - triple.t_start;
        triple.curve_3d = curve.geometry;
        // G1: a NURBS edge's end control points ARE its anchors, bit
        // for bit; T⁻¹(T(p)) lands an ulp off.
        if let (crate::geometry::curve::CurveGeometry::Nurbs(n), EdgeAnchors::Open { start, end }) =
            (&mut triple.curve_3d, &original)
        {
            if let Some(first) = n.control_points.first_mut() {
                *first = *start;
            }
            if let Some(last) = n.control_points.last_mut() {
                *last = *end;
            }
        }
        (triple.t_start, triple.t_end) = (curve.t_start, curve.t_end);
        for (pcurve, surface) in [
            (&mut triple.pcurve_a, surf_a),
            (&mut triple.pcurve_b, surf_b),
        ] {
            pcurve.shift_u(-chart_u_shift(surface, &t));
            if shift != 0.0 {
                pcurve.shift_parameter(shift);
            }
        }
        triple
    }))
}
