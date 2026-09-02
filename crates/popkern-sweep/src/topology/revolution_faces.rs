//! Topology of a FULL revolution of a closed profile (PK-408): every
//! span becomes one periodic face; adjacent faces share the circle
//! their common profile vertex sweeps, and each face closes on itself
//! along ONE seam edge (the profile span at angle 0, used by both the
//! v=0 and v=1 sides of the loop). A profile vertex on the axis sweeps
//! a point, not a circle: its circle is omitted and the two seams meet
//! at the pole vertex.
//!
//! Counts for a ring from an n-gon profile: V = n, E = 2n, F = n —
//! χ = 0, the torus the validator expects. The open-profile builder
//! used before gave every face private corners (V = E = 4n).

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use crate::error::SweepError;

/// A degree-1 pcurve in the face's (u, v) chart.
fn uv_line(u0: f64, v0: f64, u1: f64, v1: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// The v-direction boundary of a revolution surface at `u` (the circle
/// a profile point sweeps), as the surface's exact iso-curve.
fn circle_iso(surface: &NurbsSurface, u_index: usize) -> Option<NurbsCurve> {
    // Control net layout (core `revolution_surface`): for each profile
    // point, its `num_v` arc points in a row — index `u * num_v + v`.
    let mut control_points = Vec::with_capacity(surface.num_v);
    let mut weights = Vec::with_capacity(surface.num_v);
    for j in 0..surface.num_v {
        let idx = u_index * surface.num_v + j;
        control_points.push(*surface.control_points.get(idx)?);
        weights.push(*surface.weights.get(idx)?);
    }
    Some(NurbsCurve {
        degree: surface.degree_v,
        knots: surface.knots_v.clone(),
        control_points,
        weights,
    })
}

/// Distance of `point` from the axis line.
fn axis_distance(point: &Point3, axis_origin: &Point3, axis_dir: &Vec3) -> f64 {
    let d = *point - *axis_origin;
    let along = d.dot(axis_dir);
    (d - *axis_dir * along).length()
}

/// Build the welded topology of a full revolution.
///
/// `spans` are the profile spans in world coordinates with on-axis
/// spans already removed: either a closed chain (span `i` ends where
/// span `i + 1 mod n` starts; no profile vertex on the axis) or an
/// open chain whose two end vertices lie on the axis (poles).
/// `surfaces[i]` is the revolution of span `i` (u = profile, v = angle
/// normalised to [0, 1]).
pub(crate) fn build_full_revolution_topology(
    spans: &[NurbsCurve],
    surfaces: &[NurbsSurface],
    axis_origin: &Point3,
    axis_dir: &Vec3,
    tolerance: Tolerance,
) -> Result<BRep, SweepError> {
    let n = surfaces.len();
    if n == 0 || spans.len() != n {
        return Err(SweepError::Topology(
            "full revolution needs one surface per profile span".into(),
        ));
    }
    let axis_dir = axis_dir
        .normalized()
        .ok_or_else(|| SweepError::InvalidProfile("zero-length revolution axis".into()))?;
    let band = tolerance.linear() * 10.0;
    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    // Chain vertices: the start of every span, plus the end of the
    // last span when the chain does not close (it then ends on the
    // axis). A vertex on the axis sweeps itself, not a circle.
    let first_start = surfaces[0].point_at(0.0, 0.0);
    let last_end = surfaces[n - 1].point_at(1.0, 0.0);
    let closed = first_start.distance_to(&last_end) < band;
    let mut positions: Vec<Point3> = surfaces.iter().map(|s| s.point_at(0.0, 0.0)).collect();
    if !closed {
        positions.push(last_end);
    }
    let on_axis: Vec<bool> = positions
        .iter()
        .map(|p| axis_distance(p, axis_origin, &axis_dir) < band)
        .collect();
    if !closed && !(on_axis[0] && on_axis[n]) {
        return Err(SweepError::InvalidProfile(
            "a revolved profile that does not close must start and end on the axis".into(),
        ));
    }
    let vertices: Vec<_> = positions.iter().map(|p| builder.add_vertex(*p)).collect();
    let vertex_count = positions.len();

    // One closed circle edge per off-axis vertex: the v-direction
    // boundary at u = 0 of the span starting there (u = 1 of the span
    // ending there for the open chain's last vertex).
    let mut circle_fwd: Vec<Option<HalfEdgeHandle>> = vec![None; vertex_count];
    let mut circle_rev: Vec<Option<HalfEdgeHandle>> = vec![None; vertex_count];
    for j in 0..vertex_count {
        if on_axis[j] {
            continue;
        }
        let iso = if j < n {
            circle_iso(&surfaces[j], 0)
        } else {
            circle_iso(&surfaces[n - 1], surfaces[n - 1].num_u - 1)
        };
        let (edge, fwd, rev) = builder.add_edge(vertices[j], vertices[j]);
        if let Some(iso) = iso {
            builder.set_edge_curve(edge, crate::surface::canonical::line_or_nurbs_curve(&iso));
        }
        circle_fwd[j] = Some(fwd);
        circle_rev[j] = Some(rev);
    }

    for i in 0..n {
        let next = if closed { (i + 1) % n } else { i + 1 };
        // The seam: profile span i at angle 0, vertex i → vertex next.
        let (seam, seam_fwd, seam_rev) = builder.add_edge(vertices[i], vertices[next]);
        builder.set_edge_curve(
            seam,
            crate::surface::canonical::line_or_nurbs_curve(&spans[i]),
        );
        let face = builder.add_face(shell);
        builder.set_face_surface(
            face,
            SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(surfaces[i].clone()),
                orientation: true,
            },
        );
        // Loop: seam forward (v = 0, u 0→1), the next vertex's circle
        // (u = 1, v 0→1), seam backward (v = 1, u 1→0), this vertex's
        // circle backward (u = 0, v 1→0). A pole vertex contributes no
        // circle: the loop turns at the vertex.
        let mut loop_hes = vec![seam_fwd];
        builder.set_half_edge_pcurve(seam_fwd, uv_line(0.0, 0.0, 1.0, 0.0));
        if let Some(he) = circle_fwd[next] {
            loop_hes.push(he);
            builder.set_half_edge_pcurve(he, uv_line(1.0, 0.0, 1.0, 1.0));
        }
        loop_hes.push(seam_rev);
        builder.set_half_edge_pcurve(seam_rev, uv_line(1.0, 1.0, 0.0, 1.0));
        if let Some(he) = circle_rev[i] {
            loop_hes.push(he);
            builder.set_half_edge_pcurve(he, uv_line(0.0, 1.0, 0.0, 0.0));
        }
        builder
            .add_loop(&loop_hes, face)
            .map_err(|e| SweepError::Topology(format!("revolution face loop: {e}")))?;
    }
    Ok(builder.build_unchecked())
}
