//! Cap face creation for swept and lofted solids.
//!
//! Creates planar cap faces at the start and end of a sweep/loft path
//! by wiring boundary half-edges into a loop.

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{FaceHandle, HalfEdgeHandle, ShellHandle};

use crate::error::SweepError;

#[path = "caps_self_loop.rs"]
mod caps_self_loop;
use caps_self_loop::{attach_self_loop_pcurve, half_edge_runs_curve_forward};

/// Create a planar cap face from boundary half-edges.
///
/// The `half_edges` slice contains the twin half-edges of the side face
/// boundary edges that form the cap boundary. These are wired into a loop
/// and a planar surface is fitted through the boundary points.
///
/// `is_start_cap` controls the orientation of the cap normal: start caps
/// face opposite to the path direction, end caps face along it.
pub fn create_cap_face(
    builder: &mut BRepBuilder,
    shell: ShellHandle,
    half_edges: &[HalfEdgeHandle],
    cap_curves: &[NurbsCurve],
    is_start_cap: bool,
) -> Result<FaceHandle, SweepError> {
    if half_edges.is_empty() {
        return Err(SweepError::Topology("cap has no boundary edges".into()));
    }

    let face = builder.add_face(shell);

    // Compute a planar surface from the cap boundary.
    let cap_surface = compute_cap_surface(builder, half_edges, cap_curves, is_start_cap);
    builder.set_face_surface(face, cap_surface);

    builder
        .add_loop(half_edges, face)
        .map_err(|e| SweepError::Topology(format!("cap loop: {e}")))?;

    // Attach cap curve geometry to each boundary edge — with
    // canonical recognition (circular cap curves become analytic
    // Arcs; see surface::canonical).
    for (he_handle, curve) in half_edges.iter().zip(cap_curves.iter()) {
        let edge_handle = builder.brep().half_edges[*he_handle].edge;
        builder.set_edge_curve(
            edge_handle,
            crate::surface::canonical::analytic_or_nurbs_curve(curve),
        );
    }

    // Attach pcurves to cap half-edges.
    attach_cap_pcurves(builder, face, half_edges, cap_curves);

    Ok(face)
}

/// [`create_cap_face`] with inner (hole) loops — the multi-wire
/// closed sweep's caps (PK-477). The outer wire's boundary defines
/// the cap surface; each inner wire's cap-side half-edges wire an
/// inner loop with its own curves and pcurves. Inner-loop winding is
/// opposed to the outer by construction: `analyze_multi_wire`
/// normalizes inner wire profiles CW where the outer is CCW, and the
/// identical per-wire ring machinery then hands each cap chain the
/// mirrored traversal.
pub fn create_cap_face_with_inners(
    builder: &mut BRepBuilder,
    shell: ShellHandle,
    outer_half_edges: &[HalfEdgeHandle],
    outer_curves: &[NurbsCurve],
    inners: &[(Vec<HalfEdgeHandle>, Vec<NurbsCurve>)],
    is_start_cap: bool,
) -> Result<FaceHandle, SweepError> {
    let face = create_cap_face(builder, shell, outer_half_edges, outer_curves, is_start_cap)?;
    for (wire_index, (half_edges, cap_curves)) in inners.iter().enumerate() {
        if half_edges.is_empty() {
            return Err(SweepError::Topology(format!(
                "cap inner wire {wire_index} has no boundary edges"
            )));
        }
        builder
            .add_inner_loop(half_edges, face)
            .map_err(|e| SweepError::Topology(format!("cap inner loop {wire_index}: {e}")))?;
        for (he_handle, curve) in half_edges.iter().zip(cap_curves.iter()) {
            let edge_handle = builder.brep().half_edges[*he_handle].edge;
            builder.set_edge_curve(
                edge_handle,
                crate::surface::canonical::analytic_or_nurbs_curve(curve),
            );
        }
        attach_cap_pcurves(builder, face, half_edges, cap_curves);
    }
    Ok(face)
}

/// Attach pcurves to cap face half-edges.
///
/// For multi-edge caps, projects the edge curve endpoints into the
/// face's UV space and creates degree-1 line pcurves.
///
/// For self-loop caps (1 half-edge with a full-circle curve), creates
/// an exact rational circle pcurve starting at the loop vertex in the
/// half-edge's measured winding sense (see `attach_self_loop_pcurve`).
fn attach_cap_pcurves(
    builder: &mut BRepBuilder,
    face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
    cap_curves: &[NurbsCurve],
) {
    let surface = match &builder.brep().faces[face].surface {
        Some(s) => s.geometry.clone(),
        None => return,
    };

    if half_edges.len() == 1 && cap_curves.len() == 1 {
        // Self-loop cap: one full-circle edge.
        attach_self_loop_pcurve(builder, half_edges[0], &cap_curves[0], &surface);
    } else {
        // Multi-edge cap. STRAIGHT spans get exact line pcurves from
        // projected endpoints. CURVED spans (mixed-path arcs) project
        // the WHOLE control polygon: an endpoint-chord pcurve places
        // the boundary up to the sagitta inside the true face
        // (measured: a boolean's UV clip crossed a gusset cap's arc
        // pcurve at the CHORD, 7.75 mm from the true arc crossing,
        // and the intersection ring never closed). The cap chart is
        // affine, so projecting each control point with unchanged
        // weights/knots reproduces the exact rational arc in UV.
        for (&he, curve) in half_edges.iter().zip(cap_curves.iter()) {
            let (t_start, t_end) = curve.parameter_range();
            let p_start = curve.point_at(t_start);
            let p_end = curve.point_at(t_end);

            let forward = half_edge_runs_curve_forward(builder.brep(), he, curve);

            let curved = curve.degree > 1;
            if curved {
                let mut uv_cps = Vec::with_capacity(curve.control_points.len());
                let mut ok = true;
                for cp in &curve.control_points {
                    match project_point_to_uv(&surface, cp) {
                        Some((u, v)) => uv_cps.push(Point3::new(u, v, 0.0)),
                        None => {
                            ok = false;
                            break;
                        }
                    }
                }
                if ok {
                    let mut pcurve = NurbsCurve {
                        degree: curve.degree,
                        knots: curve.knots.clone(),
                        control_points: uv_cps,
                        weights: curve.weights.clone(),
                    };
                    if !forward {
                        pcurve.control_points.reverse();
                        pcurve.weights.reverse();
                        let (k0, k1) = (
                            pcurve.knots.first().copied().unwrap_or(0.0),
                            pcurve.knots.last().copied().unwrap_or(1.0),
                        );
                        pcurve.knots = pcurve.knots.iter().rev().map(|&k| k0 + k1 - k).collect();
                    }
                    builder.set_half_edge_pcurve(he, pcurve);
                    continue;
                }
            }

            let (pt_a, pt_b) = if forward {
                (p_start, p_end)
            } else {
                (p_end, p_start)
            };

            if let (Some((u0, v0)), Some((u1, v1))) = (
                project_point_to_uv(&surface, &pt_a),
                project_point_to_uv(&surface, &pt_b),
            ) {
                let pcurve = NurbsCurve {
                    degree: 1,
                    knots: vec![0.0, 0.0, 1.0, 1.0],
                    control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
                    weights: vec![1.0, 1.0],
                };
                builder.set_half_edge_pcurve(he, pcurve);
            }
        }
    }
}

/// Compute a planar surface for a cap face.
///
/// The stored plane normal is the face's OUTWARD normal, derived by
/// Newell's method over the boundary traversed in loop order. The
/// boundary is sampled from the cap CURVES (walked in each
/// half-edge's direction) so a self-loop cap — one full-circle edge,
/// whose single origin vertex makes a point-based Newell degenerate
/// — still yields a decisive loop normal. Loops wind CCW about the
/// outward normal, so no start/end flip is applied on this path; the
/// path-direction hint (`is_start_cap`) survives only as the last
/// resort when no curve geometry exists. The UV basis is built from
/// the FINAL normal so the plane frame stays right-handed
/// (`u_axis × v_axis = normal`) — a left-handed frame makes UV-CCW
/// triangulations emit inward-facing triangles.
fn compute_cap_surface(
    builder: &BRepBuilder,
    half_edges: &[HalfEdgeHandle],
    cap_curves: &[NurbsCurve],
    is_start_cap: bool,
) -> SurfaceOnFace {
    let brep = builder.brep();

    // Collect vertex positions from half-edge origins.
    let points: Vec<Point3> = half_edges
        .iter()
        .map(|&he| brep.vertices[brep.half_edges[he].origin].position)
        .collect();

    let centroid = compute_centroid(&points);

    // Primary: Newell over the curve-sampled boundary polygon in
    // loop-traversal order.
    let normal = match cap_boundary_newell(builder, half_edges, cap_curves) {
        Some(n) => n,
        None => {
            // Degenerate fallback (no curve geometry): origin-vertex
            // Newell with the historical path-direction flip.
            let n = compute_newell_normal(&points);
            if is_start_cap { n * -1.0 } else { n }
        }
    };

    // Choose u/v axes perpendicular to the FINAL normal.
    let (u_axis, v_axis) = orthonormal_basis(&normal);

    SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: centroid,
            normal,
            u_axis,
            v_axis,
        },
        orientation: true,
    }
}

/// Newell normal of the cap boundary, sampled from the cap curves and
/// traversed in loop (half-edge) order. Each curve is sampled at
/// interior resolution and reversed when the half-edge runs against
/// the curve's parameter direction. Returns `None` when the boundary
/// polygon is degenerate or curve geometry is missing.
fn cap_boundary_newell(
    builder: &BRepBuilder,
    half_edges: &[HalfEdgeHandle],
    cap_curves: &[NurbsCurve],
) -> Option<Vec3> {
    if cap_curves.len() != half_edges.len() {
        return None;
    }
    let brep = builder.brep();
    const SAMPLES: usize = 8;
    let mut polygon: Vec<Point3> = Vec::with_capacity(half_edges.len() * SAMPLES);
    for (&he, curve) in half_edges.iter().zip(cap_curves.iter()) {
        let (t0, t1) = curve.parameter_range();
        let forward = half_edge_runs_curve_forward(brep, he, curve);
        for k in 0..SAMPLES {
            let f = k as f64 / SAMPLES as f64;
            let f = if forward { f } else { 1.0 - f };
            polygon.push(curve.point_at(t0 + (t1 - t0) * f));
        }
    }
    let n = compute_newell_normal_raw(&polygon)?;
    Some(n)
}

/// Newell's method returning `None` on degenerate polygons instead of
/// a Z-up guess.
fn compute_newell_normal_raw(points: &[Point3]) -> Option<Vec3> {
    let n = points.len();
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    for i in 0..n {
        let curr = &points[i];
        let next = &points[(i + 1) % n];
        nx += (curr.y - next.y) * (curr.z + next.z);
        ny += (curr.z - next.z) * (curr.x + next.x);
        nz += (curr.x - next.x) * (curr.y + next.y);
    }
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if len < 1e-30 {
        None
    } else {
        Some(Vec3::new(nx / len, ny / len, nz / len))
    }
}

/// Compute the centroid of a set of points.
fn compute_centroid(points: &[Point3]) -> Point3 {
    let n = points.len().max(1) as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sz = 0.0;
    for p in points {
        sx += p.x;
        sy += p.y;
        sz += p.z;
    }
    Point3::new(sx / n, sy / n, sz / n)
}

/// Compute a normal using Newell's method.
fn compute_newell_normal(points: &[Point3]) -> Vec3 {
    let n = points.len();
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    for i in 0..n {
        let curr = &points[i];
        let next = &points[(i + 1) % n];
        nx += (curr.y - next.y) * (curr.z + next.z);
        ny += (curr.z - next.z) * (curr.x + next.x);
        nz += (curr.x - next.x) * (curr.y + next.y);
    }
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if len < 1e-30 {
        // Degenerate — fall back to Z-up.
        Vec3::new(0.0, 0.0, 1.0)
    } else {
        Vec3::new(nx / len, ny / len, nz / len)
    }
}

/// Build an orthonormal basis perpendicular to a given normal.
fn orthonormal_basis(normal: &Vec3) -> (Vec3, Vec3) {
    let candidate = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = normal
        .cross(&candidate)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn centroid_of_triangle() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
        ];
        let c = compute_centroid(&points);
        assert!((c.x - 1.0).abs() < 1e-10);
        assert!((c.y - 1.0).abs() < 1e-10);
        assert!(c.z.abs() < 1e-10);
    }

    #[test]
    fn newell_normal_of_xy_triangle() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let n = compute_newell_normal(&points);
        // Should point in +Z direction.
        assert!(n.z.abs() > 0.9, "expected Z-normal, got {:?}", n);
    }

    #[test]
    fn orthonormal_basis_perpendicular() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let (u, v) = orthonormal_basis(&normal);
        // u and v should be perpendicular to normal and to each other.
        assert!(normal.dot(&u).abs() < 1e-10);
        assert!(normal.dot(&v).abs() < 1e-10);
        assert!(u.dot(&v).abs() < 1e-10);
    }

    #[test]
    fn degenerate_newell_falls_back() {
        // All points are the same.
        let points = vec![Point3::new(1.0, 1.0, 1.0), Point3::new(1.0, 1.0, 1.0)];
        let n = compute_newell_normal(&points);
        // Should fall back to Z-up.
        assert!((n.z - 1.0).abs() < 1e-10);
    }
}
