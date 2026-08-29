//! T-junction boundary pre-split for face splitting (extracted from
//! `split_helpers.rs` under the file-size ratchet): split a face's
//! boundary edges at the trim curve's endpoints so `split_face` finds
//! entry/exit vertices.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Split the face's boundary edges wherever one of the trim curve's
/// endpoints lies mid-edge (no coincident vertex). Failures are
/// non-fatal: the subsequent face split reports its own precondition
/// error if the junction stays vertex-less.
pub(super) fn presplit_boundary_at_endpoints(
    brep: &mut BRep,
    face_handle: FaceHandle,
    curve_3d: &NurbsCurve,
    tolerance: &Tolerance,
) {
    use popkern_core::geometry::curve_inversion::invert_curve_at_point;
    let (t0, t1) = curve_3d.parameter_range();
    for endpoint in [curve_3d.point_at(t0), curve_3d.point_at(t1)] {
        let Some(face) = brep.faces.get(face_handle) else {
            return;
        };
        let mut loops = vec![face.outer_loop];
        loops.extend(face.inner_loops.iter().copied());
        let mut on_vertex = false;
        let mut split_target: Option<popkern_core::topology::EdgeHandle> = None;
        // PASS 1 — the COMPLETE boundary vertex scan. An existing
        // vertex within tolerance always wins over splitting an edge:
        // when a neighbouring face's trim already minted the junction
        // vertex, a single interleaved walk could reach a splittable
        // edge FIRST (loop-order dependent) and re-split at the same
        // point, minting a ~1e-13 stub edge between two coincident
        // vertices. Such a stub is a point contact, not an edge — and
        // it bridged classification flood-fill across the intersection
        // loop (measured: torus × cylinder union kept all 27
        // fragments; the y<0 junctions were saved only by loop order).
        'vertex_scan: for lh in &loops {
            let Some(loop_data) = brep.loops.get(*lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut he = start;
            loop {
                let he_data = &brep.half_edges[he];
                let vpos = brep.vertices[he_data.origin].position;
                if tolerance.points_coincident(&vpos, &endpoint) {
                    on_vertex = true;
                    break 'vertex_scan;
                }
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }
        // PASS 2 — only a boundary with NO coincident vertex may have
        // an edge split at the endpoint.
        'search: for lh in loops {
            if on_vertex {
                break;
            }
            let Some(loop_data) = brep.loops.get(lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut he = start;
            loop {
                let he_data = &brep.half_edges[he];
                let edge = he_data.edge;
                if let Some(curve) = brep.edges[edge].curve.as_ref() {
                    let inv = invert_curve_at_point(curve, endpoint, tolerance);
                    if popkern_core::trace::trace_enabled("face-split") {
                        eprintln!(
                            "    [presplit-inv] edge {edge:?} window ({:.3},{:.3}) -> {inv:?}",
                            curve.t_start, curve.t_end
                        );
                    }
                    if let Ok(t) = inv {
                        let lo = curve.t_start.min(curve.t_end);
                        let hi = curve.t_start.max(curve.t_end);
                        let span = (hi - lo).max(1e-12);
                        let margin = span * 1e-6;
                        if t > lo + margin
                            && t < hi - margin
                            && curve.geometry.point_at(t).distance_to(&endpoint)
                                < tolerance.linear() * 10.0
                        {
                            split_target = Some(edge);
                            break 'search;
                        }
                    }
                }
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }
        if !on_vertex && let Some(edge) = split_target {
            let r = brep.split_edge_at_point(edge, endpoint, tolerance);
            if popkern_core::trace::trace_enabled("face-split") {
                eprintln!(
                    "  [presplit] {face_handle:?} edge {edge:?} at ({:.1},{:.1},{:.1}) -> {r:?}",
                    endpoint.x, endpoint.y, endpoint.z
                );
            }
        } else if popkern_core::trace::trace_enabled("face-split") {
            eprintln!(
                "  [presplit] {face_handle:?} endpoint ({:.1},{:.1},{:.1}) on_vertex={on_vertex} no-target={}",
                endpoint.x,
                endpoint.y,
                endpoint.z,
                split_target.is_none()
            );
            if let Some(face) = brep.faces.get(face_handle)
                && let Some(loop_data) = brep.loops.get(face.outer_loop)
            {
                let start = loop_data.half_edge;
                let mut he = start;
                for _ in 0..12 {
                    let h = &brep.half_edges[he];
                    let v = brep.vertices[h.origin].position;
                    let w = brep.edges[h.edge].curve.as_ref().map(|c| {
                        use popkern_core::geometry::curve::CurveGeometry;
                        let kind = match &c.geometry {
                            CurveGeometry::Line { .. } => "Line".into(),
                            CurveGeometry::Circle { .. } => "Circle".into(),
                            CurveGeometry::Arc {
                                start_angle,
                                end_angle,
                                ..
                            } => format!("Arc[{start_angle:.3},{end_angle:.3}]"),
                            CurveGeometry::Nurbs(n) => format!(
                                "Nurbs[{:.3},{:.3}]",
                                n.knots.first().copied().unwrap_or(0.0),
                                n.knots.last().copied().unwrap_or(0.0)
                            ),
                            _ => "Other".into(),
                        };
                        (kind, c.t_start, c.t_end)
                    });
                    eprintln!(
                        "    [presplit-loop] he edge={:?} origin=({:.1},{:.1},{:.1}) window={w:?} d_end={:.3e}",
                        h.edge,
                        v.x,
                        v.y,
                        v.z,
                        v.distance_to(&endpoint)
                    );
                    he = h.next;
                    if he == start {
                        break;
                    }
                }
            }
        }
    }
}
