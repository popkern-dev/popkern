//! Geometric helpers for the face–face blend TrimBoth attach surgery
//! (PK-142 stage 2): trim polygons, uniform-knot polyline curves,
//! split-side witnesses, and loop utilities.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use crate::contact::ContactSolution;
use crate::error::BlendError;
use crate::face_blend_sheet::clamped_uniform_knots;

/// Build a structured attach decline.
pub(crate) fn decline(mechanism: &str, detail: String) -> BlendError {
    BlendError::AttachDeclined {
        mechanism: mechanism.into(),
        detail,
    }
}

/// The surface geometry of a face, or a structured decline.
pub(crate) fn face_surface(brep: &BRep, face: FaceHandle) -> Result<&SurfaceGeometry, BlendError> {
    brep.faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
        .ok_or_else(|| {
            decline(
                "face lookup",
                format!("face {face:?} is missing or unsurfaced"),
            )
        })
}

/// UV polygon of a face's outer loop on its own surface.
pub(crate) fn face_uv_polygon(
    brep: &BRep,
    face: FaceHandle,
) -> Result<Vec<(f64, f64)>, BlendError> {
    let surface = face_surface(brep, face)?;
    let positions = brep.face_outer_loop_positions(face);
    if positions.len() < 3 {
        return Err(decline(
            "trim polygon",
            format!("face {face:?} outer loop has {} vertices", positions.len()),
        ));
    }
    positions
        .iter()
        .map(|p| {
            surface.closest_point(p).ok_or_else(|| {
                decline(
                    "trim polygon",
                    format!("closest-point projection failed on face {face:?}"),
                )
            })
        })
        .collect()
}

/// Even-odd point-in-polygon test in UV space.
pub(crate) fn point_in_polygon(p: (f64, f64), poly: &[(f64, f64)]) -> bool {
    let mut inside = false;
    let n = poly.len();
    for i in 0..n {
        let (x1, y1) = poly[i];
        let (x2, y2) = poly[(i + 1) % n];
        if (y1 > p.1) != (y2 > p.1) {
            let x_cross = x1 + (p.1 - y1) / (y2 - y1) * (x2 - x1);
            if p.0 < x_cross {
                inside = !inside;
            }
        }
    }
    inside
}

/// Degree-1 NURBS through the points with the SAME clamped uniform
/// knot vector as the strip skin borders, so edge curves, wall
/// pcurves, and the skin all share one parameterization.
pub(crate) fn uniform_polyline_nurbs(points: &[Point3]) -> NurbsCurve {
    let n = points.len();
    let degree = 1.min(n.saturating_sub(1));
    NurbsCurve {
        degree,
        knots: clamped_uniform_knots(n, degree),
        control_points: points.to_vec(),
        weights: vec![1.0; n],
    }
}

/// UV-space variant of [`uniform_polyline_nurbs`] (z = 0).
pub(crate) fn uv_polyline_nurbs(uv: &[(f64, f64)]) -> NurbsCurve {
    let pts: Vec<Point3> = uv.iter().map(|&(u, v)| Point3::new(u, v, 0.0)).collect();
    uniform_polyline_nurbs(&pts)
}

/// Project 3D points into a surface's UV chart.
pub(crate) fn project_uv(
    surface: &SurfaceGeometry,
    points: &[Point3],
) -> Result<Vec<(f64, f64)>, BlendError> {
    points
        .iter()
        .map(|p| {
            surface.closest_point(p).ok_or_else(|| {
                decline(
                    "uv projection",
                    format!("closest-point projection failed at {p:?}"),
                )
            })
        })
        .collect()
}

/// Direction the section arc departs from one of its contact points
/// — after projection into the support tangent plane this points
/// toward the doomed (replaced) side of the wall. Mirrors the frame
/// of `section_arc_points`.
pub(crate) fn arc_departure(
    sol: &ContactSolution,
    at_first_side: bool,
    bulge: Option<Vec3>,
) -> Option<Vec3> {
    let u = sol.contact0 - sol.center;
    let to_end = sol.contact1 - sol.center;
    let normal = u.cross(&to_end);
    match normal.normalized() {
        None => {
            // Antipodal section: a half-circle bowing toward the
            // bulge hint departs both contacts along the bulge.
            if let Some(d) = bulge
                && let Some(uh) = u.normalized()
                && let Some(v) = (d - uh * d.dot(&uh)).normalized()
            {
                return Some(v);
            }
            // Degenerate straight-chord fallback.
            let chord = sol.contact1 - sol.contact0;
            let dir = if at_first_side {
                chord
            } else {
                Vec3::new(-chord.x, -chord.y, -chord.z)
            };
            dir.normalized()
        }
        Some(n) => {
            let v = n.cross(&u);
            if at_first_side {
                v.normalized()
            } else {
                let sweep = to_end.dot(&v).atan2(to_end.dot(&u));
                let sweep = if sweep >= 0.0 {
                    sweep
                } else {
                    sweep + std::f64::consts::TAU
                };
                (u * sweep.sin() - v * sweep.cos()).normalized()
            }
        }
    }
}

/// Position centroid of a face's outer-loop vertices.
pub(crate) fn loop_centroid(brep: &BRep, face: FaceHandle) -> Point3 {
    let positions = brep.face_outer_loop_positions(face);
    let n = positions.len().max(1) as f64;
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    for p in &positions {
        x += p.x;
        y += p.y;
        z += p.z;
    }
    Point3::new(x / n, y / n, z / n)
}

/// Decide which of the two split parts lies on the doomed (witness)
/// side of the split curve. Declines when the centroids do not
/// separate cleanly across the curve.
pub(crate) fn pick_doomed(
    brep: &BRep,
    part_1: FaceHandle,
    part_2: FaceHandle,
    anchor: Point3,
    witness: Vec3,
    stage: &str,
) -> Result<(FaceHandle, FaceHandle), BlendError> {
    let d1 = (loop_centroid(brep, part_1) - anchor).dot(&witness);
    let d2 = (loop_centroid(brep, part_2) - anchor).dot(&witness);
    if d1 > 0.0 && d2 < 0.0 {
        Ok((part_1, part_2))
    } else if d2 > 0.0 && d1 < 0.0 {
        Ok((part_2, part_1))
    } else {
        Err(decline(
            stage,
            format!(
                "split parts do not separate across the trim curve: \
                 witness dots {d1:.6} / {d2:.6}"
            ),
        ))
    }
}

/// Half-edges of a loop in cycle order.
pub(crate) fn loop_half_edges(brep: &BRep, lp: LoopHandle) -> Vec<HalfEdgeHandle> {
    let mut out = Vec::new();
    let start = brep.loops[lp].half_edge;
    let mut he = start;
    for _ in 0..brep.half_edges.len() + 1 {
        out.push(he);
        he = brep.half_edges[he].next;
        if he == start {
            break;
        }
    }
    out
}

/// Ensure a vertex exists at `point` on the face's outer boundary,
/// splitting the nearest boundary edge when needed.
pub(crate) fn ensure_boundary_vertex(
    work: &mut BRep,
    face: FaceHandle,
    point: Point3,
    tolerance: &Tolerance,
) -> Result<(), BlendError> {
    let lp = work.faces[face].outer_loop;
    let snap = tolerance.linear().max(1.0e-9) * 10.0;
    let hes = loop_half_edges(work, lp);
    for &he in &hes {
        let v = work.half_edges[he].origin;
        if work.vertices[v].position.distance_to(&point) < snap {
            return Ok(());
        }
    }
    // Nearest boundary edge: coarse curve sampling, then a local
    // ternary refinement (the coarse minimum alone can sit half a
    // sample spacing off the true foot — measured 5.7e-3 on a 4.1
    // span, far over the acceptance gate).
    let mut best: Option<(popkern_core::topology::EdgeHandle, f64)> = None;
    for &he in &hes {
        let edge = work.half_edges[he].edge;
        let Some(curve) = work.edges[edge].curve.as_ref() else {
            continue;
        };
        let span = curve.t_end - curve.t_start;
        let mut t_best = curve.t_start;
        let mut min_d = f64::INFINITY;
        for k in 0..=32 {
            let t = curve.t_start + span * k as f64 / 32.0;
            let d = curve.geometry.point_at(t).distance_to(&point);
            if d < min_d {
                min_d = d;
                t_best = t;
            }
        }
        let (mut lo, mut hi) = (
            (t_best - span / 32.0).max(curve.t_start.min(curve.t_end)),
            (t_best + span / 32.0).min(curve.t_start.max(curve.t_end)),
        );
        for _ in 0..40 {
            let m1 = lo + (hi - lo) / 3.0;
            let m2 = hi - (hi - lo) / 3.0;
            let d1 = curve.geometry.point_at(m1).distance_to(&point);
            let d2 = curve.geometry.point_at(m2).distance_to(&point);
            if d1 < d2 {
                hi = m2;
            } else {
                lo = m1;
            }
            min_d = min_d.min(d1.min(d2));
        }
        if best.map(|(_, bd)| min_d < bd).unwrap_or(true) {
            best = Some((edge, min_d));
        }
    }
    let gate = tolerance.linear().max(1.0e-9) * 100.0;
    let (edge, dist) = best.ok_or_else(|| {
        decline(
            "contact endpoint",
            "wall boundary has no curved edges to split".into(),
        )
    })?;
    if dist > gate {
        return Err(decline(
            "contact endpoint",
            format!("strip corner {point:?} is {dist:.3e} off the wall boundary (gate {gate:.3e})"),
        ));
    }
    work.split_edge_at_point(edge, point, tolerance)
        .map(|_| ())
        .map_err(|e| {
            decline(
                "contact endpoint",
                format!("boundary edge split failed at {point:?}: {e}"),
            )
        })
}
