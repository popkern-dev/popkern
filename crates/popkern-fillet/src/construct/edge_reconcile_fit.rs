//! Curved-fit reconciliation and vertex-snap repair helpers.
//!
//! Split from `edge_reconcile.rs` to honour the structural line
//! ratchet: `try_curved_fit` (the NURBS-involving intersection
//! approximation) and `snap_vertex_and_adjacent_chords` (the
//! vertex-move repair used by strip end arcs and corner arc
//! boundaries) are self-contained tails of that pass.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::topology::{BRep, EdgeHandle};

/// Approximate a NURBS-involving surface intersection by sampling
/// the chord between `p1` and `p2`, alternately projecting each
/// sample onto both surfaces until convergence, and fitting a
/// degree-3 NURBS through the refined points. Avoids the marching
/// solver in `intersect_surfaces` (which loops on small fillet
/// face inputs).
pub(super) fn try_curved_fit(
    brep: &mut BRep,
    edge: EdgeHandle,
    he: popkern_core::topology::HalfEdgeHandle,
    twin: popkern_core::topology::HalfEdgeHandle,
) -> bool {
    use popkern_core::geometry::fitting::fit_curve;
    use popkern_core::geometry::surface::SurfaceEval;

    let Some(he_data) = brep.half_edges.get(he) else {
        return false;
    };
    let Some(twin_data) = brep.half_edges.get(twin) else {
        return false;
    };
    let face_1 = brep.loops[he_data.loop_].face;
    let face_2 = brep.loops[twin_data.loop_].face;
    let v1 = he_data.origin;
    let v2 = twin_data.origin;
    let p1 = brep.vertices[v1].position;
    let p2 = brep.vertices[v2].position;
    let chord_len = p1.distance_to(&p2);
    if chord_len < 1e-9 {
        return false;
    }

    let surf_1 = brep.faces.get(face_1).and_then(|f| f.surface.as_ref());
    let surf_2 = brep.faces.get(face_2).and_then(|f| f.surface.as_ref());
    let (Some(sof_1), Some(sof_2)) = (surf_1, surf_2) else {
        return false;
    };

    // Sample chord, project alternately. 7 interior samples gives
    // a degree-3 fit room while staying cheap.
    let n_segments = 8;
    let mut points: Vec<Point3> = Vec::with_capacity(n_segments + 1);
    let mut params: Vec<f64> = Vec::with_capacity(n_segments + 1);
    points.push(p1);
    params.push(0.0);
    for i in 1..n_segments {
        let t = i as f64 / n_segments as f64;
        let mut p = Point3::new(
            p1.x + t * (p2.x - p1.x),
            p1.y + t * (p2.y - p1.y),
            p1.z + t * (p2.z - p1.z),
        );
        // Alternating projection converges LINEARLY near a shallow
        // (tangential) intersection — a step-size-only test with 8
        // iterations rejected exactly the strip-boundary edges the
        // pass exists for (PK-141: 2.1e-2 seam remnants on the
        // region-merged strips). Iterate deeper and accept on the
        // RESIDUAL to both surfaces, same as the swept-corner march.
        let max_iters = 200;
        let mut prev = p;
        let mut converged = false;
        for _ in 0..max_iters {
            let Some((u1, v1u)) = sof_1.geometry.closest_point(&p) else {
                return false;
            };
            p = sof_1.geometry.point_at(u1, v1u);
            let Some((u2, v2u)) = sof_2.geometry.closest_point(&p) else {
                return false;
            };
            p = sof_2.geometry.point_at(u2, v2u);
            if p.distance_to(&prev) < chord_len * 1e-5 {
                converged = true;
                break;
            }
            prev = p;
        }
        if !converged {
            // Residual acceptance: on both surfaces within the fit
            // scale even though the tangential slide has not stalled.
            let d1 = sof_1
                .geometry
                .closest_point(&p)
                .map(|(u, v)| p.distance_to(&sof_1.geometry.point_at(u, v)))
                .unwrap_or(f64::INFINITY);
            if d1 > chord_len * 1e-4 {
                return false;
            }
        }
        // Reject samples that drift far from the chord — likely
        // landed on a different intersection branch.
        let chord_pt = Point3::new(
            p1.x + t * (p2.x - p1.x),
            p1.y + t * (p2.y - p1.y),
            p1.z + t * (p2.z - p1.z),
        );
        if p.distance_to(&chord_pt) > chord_len * 0.5 {
            return false;
        }
        points.push(p);
        params.push(t);
    }
    points.push(p2);
    params.push(1.0);

    // Fit a degree-3 NURBS through the projected points.
    let degree = 3;
    let num_cps = (points.len()).clamp(degree + 1, points.len());
    let fitted = match fit_curve(&points, &params, degree, num_cps) {
        Ok(c) => c,
        Err(_) => return false,
    };

    // Verify endpoints land at p1 and p2.
    let (et_start, et_end) = fitted.parameter_range();
    let trim_start = fitted.point_at(et_start);
    let trim_end = fitted.point_at(et_end);
    let endpoint_err = trim_start.distance_to(&p1).max(trim_end.distance_to(&p2));
    if endpoint_err > (chord_len * 0.05).max(1e-6) {
        return false;
    }

    if let Some(edge_mut) = brep.edges.get_mut(edge) {
        edge_mut.curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(fitted),
            t_start: et_start,
            t_end: et_end,
        });
        return true;
    }
    false
}

#[allow(dead_code)]
fn pt_zero() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}

/// Move `vertex` to `target` and rebuild adjacent edges' endpoint
/// geometry: chord Lines are rebuilt between their (moved) vertex
/// positions; NURBS/polyline curves get their nearest endpoint
/// control point snapped when it sits within the same class of
/// distance. `except` (the edge being re-trimmed) is left alone.
pub(crate) fn snap_vertex_and_adjacent_chords(
    brep: &mut BRep,
    vertex: popkern_core::topology::VertexHandle,
    target: popkern_core::geometry::Point3,
    except: EdgeHandle,
) {
    use popkern_core::geometry::curve::CurveEval;
    let Some(v) = brep.vertices.get_mut(vertex) else {
        return;
    };
    let old_pos = v.position;
    if old_pos.distance_to(&target) < 1e-12 {
        return;
    }
    v.position = target;

    // Collect edges incident to the vertex.
    let incident: Vec<EdgeHandle> = brep
        .edges
        .iter()
        .filter(|(eh, e)| {
            *eh != except && {
                let he = e.half_edge;
                brep.half_edges.get(he).is_some_and(|h| {
                    h.origin == vertex
                        || brep
                            .half_edges
                            .get(h.twin)
                            .is_some_and(|t| t.origin == vertex)
                })
            }
        })
        .map(|(eh, _)| eh)
        .collect();
    for eh in incident {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let Some(curve) = edge.curve.clone() else {
            continue;
        };
        let he = edge.half_edge;
        let (va, vb) = {
            let h = &brep.half_edges[he];
            (h.origin, brep.half_edges[h.twin].origin)
        };
        let pa = brep.vertices[va].position;
        let pb = brep.vertices[vb].position;
        match &curve.geometry {
            popkern_core::geometry::CurveGeometry::Line { .. } => {
                // Rebuild the line between the current vertex
                // positions.
                let dir = pb - pa;
                let len = dir.length();
                if len > 1e-12
                    && let Some(e) = brep.edges.get_mut(eh)
                {
                    e.curve = Some(popkern_core::geometry::CurveOnEdge {
                        geometry: popkern_core::geometry::CurveGeometry::Line {
                            origin: pa,
                            direction: dir * (1.0 / len),
                        },
                        t_start: 0.0,
                        t_end: len,
                    });
                }
            }
            popkern_core::geometry::CurveGeometry::Nurbs(n) => {
                // Snap the endpoint control point nearest the moved
                // vertex.
                let ceval: &dyn CurveEval = &curve.geometry;
                let start = ceval.point_at(curve.t_start);
                let end = ceval.point_at(curve.t_end);
                let mut n2 = n.clone();
                let snap_first = start.distance_to(&old_pos) <= end.distance_to(&old_pos);
                let moved = if snap_first {
                    n2.control_points.first_mut()
                } else {
                    n2.control_points.last_mut()
                };
                if let Some(cp) = moved
                    && cp.distance_to(&old_pos) < 1e-6
                {
                    *cp = target;
                    if let Some(e) = brep.edges.get_mut(eh) {
                        e.curve = Some(popkern_core::geometry::CurveOnEdge {
                            geometry: popkern_core::geometry::CurveGeometry::Nurbs(n2),
                            t_start: curve.t_start,
                            t_end: curve.t_end,
                        });
                    }
                }
            }
            _ => {
                // Analytic curves (arcs etc.) keep their geometry;
                // a moved vertex within the snap class stays within
                // downstream association gates.
            }
        }
    }
}
