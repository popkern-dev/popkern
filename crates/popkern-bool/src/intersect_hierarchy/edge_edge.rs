//! Edge-edge intersection detection between two BRep bodies.
//!
//! Two edges intersect when their underlying curves cross at a point
//! (or overlap along a segment — handled separately by coincidence detection).
//! Uses `intersect_curves` from popkern-core for the actual computation.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::{CurveEval, CurveIntersection, intersect_curves};
use popkern_core::topology::BRep;

use super::result::EdgeEdgeHit;

/// Find all edge-edge intersection points between bodies A and B.
pub(crate) fn detect_edge_edge(
    a: &BRep,
    b: &BRep,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Vec<EdgeEdgeHit> {
    let mut hits = Vec::new();

    for (ea_h, ea) in &a.edges {
        let ca = match &ea.curve {
            Some(c) => c,
            None => continue,
        };
        let a_range = (ca.t_start.min(ca.t_end), ca.t_start.max(ca.t_end));
        let aabb_a = match edge_aabb(a, ea_h) {
            Some(bb) => bb,
            None => continue,
        };

        for (eb_h, eb) in &b.edges {
            let cb = match &eb.curve {
                Some(c) => c,
                None => continue,
            };
            let b_range = (cb.t_start.min(cb.t_end), cb.t_start.max(cb.t_end));
            let aabb_b = match edge_aabb(b, eb_h) {
                Some(bb) => bb,
                None => continue,
            };

            // AABB filter
            if !aabb_a.overlaps(&aabb_b, tolerance.linear) {
                continue;
            }

            // Intersect the two curves
            let result = intersect_curves(&ca.geometry, a_range, &cb.geometry, b_range, tolerance);
            if let Ok(result) = result {
                for isect in &result.intersections {
                    if let CurveIntersection::Point { t1, t2, point, .. } = isect {
                        // Filter out hits at endpoints (those are VV or VE)
                        let a_len = (a_range.1 - a_range.0).abs();
                        let b_len = (b_range.1 - b_range.0).abs();
                        let eps_a = a_len * 0.01;
                        let eps_b = b_len * 0.01;
                        let at_a_end =
                            (*t1 - a_range.0).abs() < eps_a || (*t1 - a_range.1).abs() < eps_a;
                        let at_b_end =
                            (*t2 - b_range.0).abs() < eps_b || (*t2 - b_range.1).abs() < eps_b;
                        if at_a_end || at_b_end {
                            continue;
                        }

                        hits.push(EdgeEdgeHit {
                            edge_a: ea_h,
                            edge_b: eb_h,
                            t_a: *t1,
                            t_b: *t2,
                            point: *point,
                        });
                    }
                }
            }
        }
    }

    hits
}

/// Compute AABB for an edge from its endpoint vertex positions.
fn edge_aabb(brep: &BRep, edge_handle: popkern_core::topology::EdgeHandle) -> Option<Aabb> {
    let edge = &brep.edges[edge_handle];
    let he = brep.half_edges.get(edge.half_edge)?;
    let v_start = brep.vertices.get(he.origin)?.position;
    let twin = brep.half_edges.get(he.twin)?;
    let v_end = brep.vertices.get(twin.origin)?.position;

    // For curved edges, also sample the midpoint
    if let Some(curve) = &edge.curve {
        let t_mid = (curve.t_start + curve.t_end) * 0.5;
        let p_mid = curve.geometry.point_at(t_mid);
        Aabb::from_points(&[v_start, v_end, p_mid])
    } else {
        Aabb::from_points(&[v_start, v_end])
    }
}
