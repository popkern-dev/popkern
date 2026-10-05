//! Vertex snap onto recomputed analytic edge curves (draft stage
//! 3.5): after edge-curve recomputation the new curves are the
//! geometric authority; endpoint vertices and parameter ranges are
//! reconciled onto them here.

use popkern_math::PortableMath;
use std::collections::HashMap;

use crate::error::DraftWarning;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// Snap the endpoint vertices of recomputed edges onto their new
/// curves, then re-derive each recomputed edge's parameter range from
/// the snapped positions. Vertex recomputation only solves plane
/// corners, so rim vertices of a drafted cylinder→cone still sit at
/// the OLD radius when the edges are recomputed; the new analytic
/// curves are the authority. A vertex shared by several recomputed
/// curves keeps the first snap and records a warning if another curve
/// disagrees beyond 10× linear tolerance.
pub fn snap_vertices_to_recomputed_curves(
    brep: &mut BRep,
    recomputed: &[EdgeHandle],
    tolerance: &Tolerance,
    warnings: &mut Vec<DraftWarning>,
) {
    use popkern_core::topology::VertexHandle;
    let mut snapped: HashMap<VertexHandle, popkern_core::geometry::Point3> = HashMap::new();
    // Circles first: a rim circle lies in its cap plane, so the exact
    // corner point IS on the circle — a line processed first would
    // instead park the vertex at its perpendicular foot, off the cap
    // by O(r·tan²α) (measured: top seam vertex z 115.963 vs 116).
    let mut ordered: Vec<EdgeHandle> = recomputed.to_vec();
    ordered.sort_by_key(|&eh| {
        let is_circle = brep
            .edges
            .get(eh)
            .and_then(|e| e.curve.as_ref())
            .is_some_and(|c| matches!(c.geometry, CurveGeometry::Circle { .. }));
        if is_circle { 0 } else { 1 }
    });
    for &eh in &ordered {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let Some(curve) = edge.curve.clone() else {
            continue;
        };
        let he = &brep.half_edges[edge.half_edge];
        let twin_origin = brep.half_edges[he.twin].origin;
        for vh in [he.origin, twin_origin] {
            let pos = brep.vertices[vh].position;
            let target = closest_point_on_curve(&curve.geometry, pos);
            match snapped.get(&vh) {
                Some(prev) => {
                    let dev = prev.distance_to(&target);
                    if dev > 10.0 * tolerance.linear() {
                        warnings.push(DraftWarning::VertexSnapDisagreement {
                            vertex: vh,
                            deviation: dev,
                        });
                    }
                }
                None => {
                    snapped.insert(vh, target);
                    brep.vertices[vh].position = target;
                }
            }
        }
    }
    // Re-derive parameter ranges from the snapped endpoints.
    for &eh in recomputed {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let he_handle = edge.half_edge;
        let he = &brep.half_edges[he_handle];
        let origin = brep.vertices[he.origin].position;
        let twin_origin_v = brep.half_edges[he.twin].origin;
        let closed = he.origin == twin_origin_v;
        let end = brep.vertices[twin_origin_v].position;
        let Some(curve) = brep.edges[eh].curve.clone() else {
            continue;
        };
        let (t_start, t_end) = match &curve.geometry {
            CurveGeometry::Line {
                origin: lo,
                direction,
            } => {
                let len2 = direction.dot(direction);
                if len2 < tolerance.numerical_floor() {
                    continue;
                }
                (
                    (origin - *lo).dot(direction) / len2,
                    (end - *lo).dot(direction) / len2,
                )
            }
            CurveGeometry::Circle { center, normal, .. } => {
                let (p1, p2) = popkern_core::geometry::analytic_curves::circle_basis(*normal);
                let d = origin - *center;
                let a0 = d.dot(&p2).patan2(d.dot(&p1));
                if closed {
                    (a0, a0 + std::f64::consts::TAU)
                } else {
                    let de = end - *center;
                    let mut a1 = de.dot(&p2).patan2(de.dot(&p1));
                    while a1 <= a0 {
                        a1 += std::f64::consts::TAU;
                    }
                    (a0, a1)
                }
            }
            _ => continue,
        };
        if let Some(edge) = brep.edges.get_mut(eh)
            && let Some(c) = edge.curve.as_mut()
        {
            c.t_start = t_start;
            c.t_end = t_end;
        }
    }
}

/// Closest point on an analytic curve (Line or Circle) to a query
/// position. Other variants return the query unchanged.
fn closest_point_on_curve(
    geometry: &CurveGeometry,
    pos: popkern_core::geometry::Point3,
) -> popkern_core::geometry::Point3 {
    match geometry {
        CurveGeometry::Line { origin, direction } => {
            let len2 = direction.dot(direction);
            if len2 <= 0.0 {
                return pos;
            }
            *origin + *direction * ((pos - *origin).dot(direction) / len2)
        }
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => {
            let d = pos - *center;
            let radial = d - *normal * d.dot(normal);
            match radial.normalized() {
                Some(dir) => *center + dir * *radius,
                None => pos,
            }
        }
        _ => pos,
    }
}
