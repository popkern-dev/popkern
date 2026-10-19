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
//!
//! Geometry (PK-409): when EVERY span is recognised as a canonical
//! surface of revolution (plane, cylinder, cone, torus, sphere), the
//! faces carry the exact surface, pcurves are authored in its native
//! chart, and circle / seam edges are angle-parameterised `Arc`s —
//! co-parameterised with the angle-linear pcurves, so the exact
//! revolution meshers serve seams and poles. Otherwise every face
//! stays NURBS with iso-curve edges (co-parameterised with the
//! rational surface parameter). Mixing the two conventions on one
//! shared edge is deliberately not attempted.

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::geometry::analytic_curves::circle_basis;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle};
use popkern_math::PortableMath;
use std::f64::consts::TAU;

use crate::error::SweepError;
use crate::surface::canonical::canonical_normal_agreement;
use crate::topology::revolution_canonical::{
    RevolutionChart, RevolvedSurface, recognize_revolved_span,
};

/// A degree-1 pcurve between two chart points.
fn uv_line(u0: f64, v0: f64, u1: f64, v1: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// A rational full circle in a planar chart, starting at azimuth
/// `phi0` and running in `sense` (four quadrant arcs).
fn planar_circle(center: (f64, f64), radius: f64, phi0: f64, sense: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let at = |angle: f64, scale: f64| {
        Point3::new(
            center.0 + scale * radius * angle.pcos(),
            center.1 + scale * radius * angle.psin(),
            0.0,
        )
    };
    let mut control_points = Vec::with_capacity(9);
    let mut weights = Vec::with_capacity(9);
    for k in 0..4 {
        let a0 = phi0 + sense * (k as f64) * std::f64::consts::FRAC_PI_2;
        let a1 = a0 + sense * std::f64::consts::FRAC_PI_4;
        control_points.push(at(a0, 1.0));
        weights.push(1.0);
        control_points.push(at(a1, std::f64::consts::SQRT_2));
        weights.push(w);
    }
    control_points.push(at(phi0, 1.0));
    weights.push(1.0);
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points,
        weights,
    }
}

/// The v-direction boundary of a revolution surface at `u_index`.
fn circle_iso(surface: &NurbsSurface, u_index: usize) -> Option<NurbsCurve> {
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

fn axis_distance(point: &Point3, axis_origin: &Point3, axis_dir: &Vec3) -> f64 {
    let d = *point - *axis_origin;
    let along = d.dot(axis_dir);
    (d - *axis_dir * along).length()
}

/// The exact circle a vertex sweeps: an `Arc` anchored at the vertex,
/// running in the sweep's sense about the axis.
fn vertex_circle(point: &Point3, axis_origin: &Point3, axis_dir: &Vec3, sense: f64) -> CurveOnEdge {
    let d = *point - *axis_origin;
    let along = d.dot(axis_dir);
    let center = *axis_origin + *axis_dir * along;
    let normal = *axis_dir * sense;
    let radius = (d - *axis_dir * along).length();
    let (e1, e2) = circle_basis(normal);
    let r = *point - center;
    let theta0 = r.dot(&e2).patan2(r.dot(&e1));
    CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle: theta0,
            end_angle: theta0 + TAU,
        },
        t_start: theta0,
        t_end: theta0 + TAU,
    }
}

/// The sense of the sweep about the axis (`+1` = right-handed).
fn sweep_sense(surface: &NurbsSurface, axis_origin: &Point3, axis_dir: &Vec3) -> f64 {
    let p0 = surface.point_at(0.5, 0.0);
    let p1 = surface.point_at(0.5, 0.05);
    let r0 = p0 - *axis_origin;
    let r0 = r0 - *axis_dir * r0.dot(axis_dir);
    let r1 = p1 - *axis_origin;
    let r1 = r1 - *axis_dir * r1.dot(axis_dir);
    if r0.cross(&r1).dot(axis_dir) >= 0.0 {
        1.0
    } else {
        -1.0
    }
}

/// Pcurves of one canonical face's loop: `(seam_fwd, circle_next,
/// seam_rev, circle_this)`.
fn canonical_pcurves(chart: &RevolutionChart) -> (NurbsCurve, NurbsCurve, NurbsCurve, NurbsCurve) {
    match chart {
        RevolutionChart::Angular {
            phi0,
            sense,
            sweep,
            v_a,
            v_b,
        } => {
            let phi1 = phi0 + sense * sweep;
            (
                uv_line(*phi0, *v_a, *phi0, *v_b),
                uv_line(*phi0, *v_b, phi1, *v_b),
                uv_line(phi1, *v_b, phi1, *v_a),
                uv_line(phi1, *v_a, *phi0, *v_a),
            )
        }
        RevolutionChart::Planar {
            center,
            phi0,
            sense,
            r_a,
            r_b,
        } => {
            let at = |r: f64| (center.0 + r * phi0.pcos(), center.1 + r * phi0.psin());
            let (a, b) = (at(*r_a), at(*r_b));
            (
                uv_line(a.0, a.1, b.0, b.1),
                planar_circle(*center, *r_b, *phi0, *sense),
                uv_line(b.0, b.1, a.0, a.1),
                planar_circle(*center, *r_a, *phi0, -sense),
            )
        }
    }
}

/// The face chart of a full revolution's construction surface: its
/// transpose (angle in u), in periodic form in u when the angle rows
/// close within the linear tolerance; otherwise the clamped transpose.
fn periodic_angle_chart(surface: &NurbsSurface, tolerance: Tolerance) -> NurbsSurface {
    let chart = surface.transposed();
    chart
        .periodic_u_from_closed_clamped(tolerance.linear())
        .unwrap_or(chart)
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
    // Canonical geometry only when every span is recognised.
    let recognized: Vec<Option<RevolvedSurface>> = spans
        .iter()
        .zip(surfaces)
        .map(|(span, surface)| {
            recognize_revolved_span(
                span,
                surface,
                axis_origin,
                &axis_dir,
                TAU,
                tolerance.linear(),
            )
        })
        .collect();
    let canonical = recognized.iter().all(Option::is_some);
    let sense = sweep_sense(&surfaces[0], axis_origin, &axis_dir);

    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

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
    // Planar faces (annuli, discs) carry no seam slit: a plane with a
    // seam edge on both of its sides is an open sheet to the boolean
    // and a slit to every mesher. Their loops are circles only, so a
    // pole vertex exists only where a non-planar span's seam ends.
    let planar: Vec<bool> = recognized
        .iter()
        .map(|r| {
            canonical
                && matches!(
                    r,
                    Some(RevolvedSurface {
                        chart: RevolutionChart::Planar { .. },
                        ..
                    })
                )
        })
        .collect();
    let vertex_count = positions.len();
    let seam_touches = |j: usize| -> bool {
        let starts = j < n && !planar[j];
        let ends = if closed {
            !planar[(j + n - 1) % n]
        } else {
            j > 0 && !planar[j - 1]
        };
        starts || ends
    };
    let vertices: Vec<Option<_>> = positions
        .iter()
        .enumerate()
        .map(|(j, p)| (!on_axis[j] || seam_touches(j)).then(|| builder.add_vertex(*p)))
        .collect();

    let mut circle_fwd: Vec<Option<HalfEdgeHandle>> = vec![None; vertex_count];
    let mut circle_rev: Vec<Option<HalfEdgeHandle>> = vec![None; vertex_count];
    for j in 0..vertex_count {
        let (false, Some(vertex)) = (on_axis[j], vertices[j]) else {
            continue;
        };
        let (edge, fwd, rev) = builder.add_edge(vertex, vertex);
        if canonical {
            builder.set_edge_curve(
                edge,
                vertex_circle(&positions[j], axis_origin, &axis_dir, sense),
            );
        } else {
            let iso = if j < n {
                circle_iso(&surfaces[j], 0)
            } else {
                circle_iso(&surfaces[n - 1], surfaces[n - 1].num_u - 1)
            };
            if let Some(iso) = iso {
                builder.set_edge_curve(edge, crate::surface::canonical::line_or_nurbs_curve(&iso));
            }
        }
        circle_fwd[j] = Some(fwd);
        circle_rev[j] = Some(rev);
    }

    for i in 0..n {
        let next = if closed { (i + 1) % n } else { i + 1 };
        let face = builder.add_face(shell);
        let (surface_on_face, pcurves) = match recognized[i].as_ref().filter(|_| canonical) {
            Some(exact) => (
                SurfaceOnFace {
                    orientation: canonical_normal_agreement(&surfaces[i], &exact.geometry),
                    geometry: exact.geometry.clone(),
                },
                canonical_pcurves(&exact.chart),
            ),
            None => (
                // The face chart is the TRANSPOSE of the construction
                // surface: angle in u, profile in v (PK-456). The
                // closed direction then matches every other closed
                // chart in the kernel — analytic revolutions, swept
                // and lofted walls — so the boolean chart walk and
                // the closed-chart meshers apply instead of the
                // v-closed decline (PK-455). Transposition reverses
                // the chart normal; the flipped orientation keeps the
                // material side. The angle direction is the exact
                // rational circle, smooth across its seam, so the wall
                // is written in periodic form there: a ring with no
                // seam boundary in its chart (PK-893).
                SurfaceOnFace {
                    geometry: SurfaceGeometry::Nurbs(periodic_angle_chart(&surfaces[i], tolerance)),
                    orientation: false,
                },
                (
                    uv_line(0.0, 0.0, 0.0, 1.0),
                    uv_line(0.0, 1.0, 1.0, 1.0),
                    uv_line(1.0, 1.0, 1.0, 0.0),
                    uv_line(1.0, 0.0, 0.0, 0.0),
                ),
            ),
        };
        builder.set_face_surface(face, surface_on_face);
        if planar[i] {
            // Outer loop: the larger circle; the other is the hole.
            let r_this = axis_distance(&positions[i], axis_origin, &axis_dir);
            let r_next = axis_distance(&positions[next], axis_origin, &axis_dir);
            let this_loop = circle_rev[i].map(|he| (he, pcurves.3.clone()));
            let next_loop = circle_fwd[next].map(|he| (he, pcurves.1.clone()));
            let (outer, inner) = if r_next >= r_this {
                (next_loop, this_loop)
            } else {
                (this_loop, next_loop)
            };
            let Some((outer_he, outer_pc)) = outer else {
                return Err(SweepError::Topology(
                    "planar revolution face without a circle boundary".into(),
                ));
            };
            builder.set_half_edge_pcurve(outer_he, outer_pc);
            builder
                .add_loop(&[outer_he], face)
                .map_err(|e| SweepError::Topology(format!("planar revolution face loop: {e}")))?;
            if let Some((inner_he, inner_pc)) = inner {
                builder.set_half_edge_pcurve(inner_he, inner_pc);
                builder
                    .add_inner_loop(&[inner_he], face)
                    .map_err(|e| SweepError::Topology(format!("planar revolution hole: {e}")))?;
            }
            continue;
        }
        let (Some(v_i), Some(v_next)) = (vertices[i], vertices[next]) else {
            return Err(SweepError::Topology(
                "revolution seam without its vertices".into(),
            ));
        };
        let (seam, seam_fwd, seam_rev) = builder.add_edge(v_i, v_next);
        builder.set_edge_curve(
            seam,
            if canonical {
                crate::surface::canonical::analytic_or_nurbs_curve(&spans[i])
            } else {
                crate::surface::canonical::line_or_nurbs_curve(&spans[i])
            },
        );
        let mut loop_hes = vec![seam_fwd];
        builder.set_half_edge_pcurve(seam_fwd, pcurves.0);
        if let Some(he) = circle_fwd[next] {
            loop_hes.push(he);
            builder.set_half_edge_pcurve(he, pcurves.1);
        }
        loop_hes.push(seam_rev);
        builder.set_half_edge_pcurve(seam_rev, pcurves.2);
        if let Some(he) = circle_rev[i] {
            loop_hes.push(he);
            builder.set_half_edge_pcurve(he, pcurves.3);
        }
        builder
            .add_loop(&loop_hes, face)
            .map_err(|e| SweepError::Topology(format!("revolution face loop: {e}")))?;
    }
    Ok(builder.build_unchecked())
}
