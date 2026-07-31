//! Ruled fan cap — the always-constructible hole filler.
//!
//! Per hole side, the sector surface is the ruling between the pole
//! and the side's own rim curve: `S(u, v) = (1−v)·pole + v·rim(u)`,
//! with `u` in the rim curve's own parameter. Every boundary curve of
//! the committed fan face lies on its sector **by construction** —
//! the rim is the `v = 1` row verbatim, and the spokes are the
//! straight `u = const` rulings shared exactly with the neighbouring
//! sector — so all pcurves are exact parameter lines and the commit
//! needs no projection at all.
//!
//! Continuity is G⁰ across the spokes (a fan of rulings creases like
//! a cone net). That is the honest healing baseline: it closes any
//! detected hole whose edges carry curves, for any rim shape, and the
//! [`super::hole_heal`] report records when it served because the G¹
//! fillers declined.

use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{CurveOnEdge, NurbsSurface, Point3, SurfaceGeometry};
use crate::topology::BRep;

use super::hole_commit::HoleCommitError;
use super::hole_commit_geometry::{HoleCapChart, HoleCapGeometry, HoleCapSectorGeometry};
use super::hole_fill::DetectedHole;

/// Build a ruled fan cap for a detected hole.
///
/// The pole is the centroid of the hole's ring vertices. Each sector
/// rules between the pole and one rim edge's curve, taken in the
/// edge's own parameterization so the committed rim pcurve is exactly
/// linear.
///
/// # Errors
///
/// [`HoleCommitError::UnsupportedSurfaceForm`] when a rim edge has no
/// curve to rule against or a spoke ruling is degenerate (pole
/// coincides with a ring vertex within `tolerance`).
pub fn ruled_fan_cap(
    brep: &BRep,
    hole: &DetectedHole,
    tolerance: f64,
) -> Result<HoleCapGeometry, HoleCommitError> {
    let mut rims = Vec::with_capacity(hole.boundary_edges.len());
    for (index, &edge) in hole.boundary_edges.iter().enumerate() {
        let curve = brep
            .edges
            .get(edge)
            .and_then(|edge_data| edge_data.curve.as_ref())
            .ok_or_else(|| HoleCommitError::UnsupportedSurfaceForm {
                sector: index,
                detail: "rim edge has no curve to rule against".into(),
            })?;
        let rim = curve
            .geometry
            .to_nurbs_in_range((curve.t_start, curve.t_end));
        let start = rim.control_points.first().copied().ok_or_else(|| {
            HoleCommitError::UnsupportedSurfaceForm {
                sector: index,
                detail: "rim curve has no control points".into(),
            }
        })?;
        let end = rim.control_points.last().copied().ok_or_else(|| {
            HoleCommitError::UnsupportedSurfaceForm {
                sector: index,
                detail: "rim curve has no control points".into(),
            }
        })?;
        rims.push((rim, start, end));
    }

    // Pole: centroid of the distinct ring vertex positions. Rim curve
    // endpoints are the ring vertices (each appears twice — once as an
    // end, once as the next start), so averaging the starts covers
    // every vertex exactly once whichever way each edge is oriented.
    let inverse = 1.0 / rims.len() as f64;
    let mut pole = Point3::new(0.0, 0.0, 0.0);
    for (_, start, end) in &rims {
        let mid = Point3::new(
            0.5 * (start.x + end.x),
            0.5 * (start.y + end.y),
            0.5 * (start.z + end.z),
        );
        // Side midpoints of endpoints reproduce the vertex centroid on
        // a closed ring (each vertex contributes half twice).
        pole = Point3::new(
            pole.x + mid.x * inverse,
            pole.y + mid.y * inverse,
            pole.z + mid.z * inverse,
        );
    }

    let mut sectors = Vec::with_capacity(rims.len());
    for (index, (rim, start, end)) in rims.into_iter().enumerate() {
        let spoke = |target: Point3| {
            CurveOnEdge::line_segment(pole, target, tolerance).ok_or_else(|| {
                HoleCommitError::UnsupportedSurfaceForm {
                    sector: index,
                    detail: format!(
                        "spoke ruling from pole {pole:?} to ring vertex {target:?} is degenerate"
                    ),
                }
            })
        };
        let spoke_at_rim_start = spoke(start)?;
        let spoke_at_rim_end = spoke(end)?;
        sectors.push(HoleCapSectorGeometry {
            surface: SurfaceGeometry::Nurbs(ruled_to_pole(&rim, pole)),
            rim_start: start,
            rim_end: end,
            spoke_at_rim_start,
            spoke_at_rim_end,
        });
    }
    Ok(HoleCapGeometry {
        pole,
        sectors,
        chart: HoleCapChart::RuledFanAnalytic,
    })
}

/// Ruling between a point and a NURBS curve: `v = 0` row collapsed at
/// the pole (with the rim's weights, so rational rims rule exactly),
/// `v = 1` row the rim verbatim, `u` in the rim's own parameter.
fn ruled_to_pole(rim: &crate::geometry::NurbsCurve, pole: Point3) -> NurbsSurface {
    let count = rim.control_points.len();
    let mut control_points = Vec::with_capacity(count * 2);
    let mut weights = Vec::with_capacity(count * 2);
    for (point, weight) in rim.control_points.iter().zip(&rim.weights) {
        control_points.push(pole);
        control_points.push(*point);
        weights.push(*weight);
        weights.push(*weight);
    }
    NurbsSurface {
        degree_u: rim.degree,
        degree_v: 1,
        knots_u: rim.knots.clone(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: count,
        num_v: 2,
        weights,
    }
}

/// Author the exact parameter-line pcurves for one committed ruled
/// sector face: rim on the `v = 1` row at the edge's own parameters,
/// spokes on the `u = const` rulings.
pub(crate) fn author_ruled_sector_pcurves(
    brep: &mut BRep,
    face: crate::topology::FaceHandle,
    pole_vertex: crate::topology::VertexHandle,
    tolerance: f64,
) -> Result<(), HoleCommitError> {
    let Some(surface) = brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|on_face| on_face.geometry.clone())
    else {
        return Ok(());
    };
    let ((u_min, u_max), _) = surface.parameter_range();
    let Some(loop_data) = brep
        .faces
        .get(face)
        .and_then(|face_data| brep.loops.get(face_data.outer_loop))
    else {
        return Ok(());
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    loop {
        let Some((next, origin, edge)) = brep
            .half_edges
            .get(cursor)
            .map(|half_edge| (half_edge.next, half_edge.origin, half_edge.edge))
        else {
            return Ok(());
        };
        let destination = brep.half_edges.get(next).map(|half_edge| half_edge.origin);
        let is_spoke = origin == pole_vertex || destination == Some(pole_vertex);
        let curve_range = brep
            .edges
            .get(edge)
            .and_then(|edge_data| edge_data.curve.as_ref())
            .map(|curve| (curve.t_start, curve.t_end));
        if let Some((t_start, t_end)) = curve_range {
            let pcurve = if is_spoke {
                // Which ruling: the junction vertex sits at one end of
                // the rim row.
                let junction_vertex = if origin == pole_vertex {
                    destination
                } else {
                    Some(origin)
                };
                let junction = junction_vertex
                    .and_then(|vertex| brep.vertices.get(vertex))
                    .map(|vertex| vertex.position);
                let Some(junction) = junction else {
                    cursor = next;
                    if cursor == start {
                        break;
                    }
                    continue;
                };
                let at_u_min = surface.point_at(u_min, 1.0).distance_to(&junction);
                let at_u_max = surface.point_at(u_max, 1.0).distance_to(&junction);
                let ruling_u = if at_u_min <= at_u_max { u_min } else { u_max };
                if at_u_min.min(at_u_max) > tolerance {
                    return Err(HoleCommitError::EndpointMismatch {
                        sector: 0,
                        role: "spoke ruling junction",
                        distance: at_u_min.min(at_u_max),
                        budget: tolerance,
                    });
                }
                // Spoke curves run pole (t_start) → junction (t_end),
                // and the ruling runs pole (v=0) → rim (v=1).
                Pcurve::linear((ruling_u, 0.0), (ruling_u, 1.0), t_start, t_end)
            } else {
                // Rim: the surface's u axis IS the edge curve's own
                // parameter (the ruling was built from this curve).
                Pcurve::linear((t_start, 1.0), (t_end, 1.0), t_start, t_end)
            };
            brep.half_edges[cursor].pcurve = Some(pcurve);
        }
        cursor = next;
        if cursor == start {
            break;
        }
    }
    Ok(())
}
