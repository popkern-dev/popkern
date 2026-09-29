//! Queryable facts about a body's edges and faces (split from
//! `select.rs` under the file-size ratchet).

use super::{CurveKind, EdgeInfo, FaceInfo, SurfaceKind};
use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;

/// Collect the queryable facts about every edge of a body.
pub fn edge_infos(brep: &BRep) -> Vec<EdgeInfo> {
    let mut out = Vec::new();
    for (eh, edge) in brep.edges.iter() {
        let (curve, midpoint, direction, chord_length) = match edge.curve.as_ref() {
            Some(c) => {
                let ceval: &dyn CurveEval = &c.geometry;
                let mid = ceval.point_at(0.5 * (c.t_start + c.t_end));
                let p0 = ceval.point_at(c.t_start);
                let p1 = ceval.point_at(c.t_end);
                let chord = p1 - p0;
                let len = chord.length();
                let dir = if len > 1e-9 {
                    Some(chord * (1.0 / len))
                } else {
                    None
                };
                let kind = match &c.geometry {
                    CurveGeometry::Line { .. } => CurveKind::Line,
                    CurveGeometry::Circle { .. } | CurveGeometry::Arc { .. } => CurveKind::Circular,
                    CurveGeometry::Ellipse { .. }
                    | CurveGeometry::Parabola { .. }
                    | CurveGeometry::Hyperbola { .. } => CurveKind::Conic,
                    // Boolean output edges carry their intersection
                    // records as NURBS even when the geometry is a
                    // straight line (a plane×plane cut). Classify by
                    // GEOMETRY, not representation: a curve whose
                    // interior samples lie on the endpoint chord is a
                    // line — otherwise `straight()` silently skips
                    // exactly the edges a dress-up pass wants
                    // (measured: an arm's rest-face side edges were
                    // unselectable while the same cut's vertical
                    // edges, stored analytically, matched).
                    CurveGeometry::Nurbs(_) => {
                        let straight = len > 1e-9 && {
                            let dirv = chord * (1.0 / len);
                            (1..8).all(|i| {
                                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                                let p = ceval.point_at(t);
                                let v = p - p0;
                                let along = v.dot(&dirv);
                                (v - dirv * along).length() < 1e-9 * (1.0 + len)
                            })
                        };
                        if straight {
                            CurveKind::Line
                        } else {
                            CurveKind::Nurbs
                        }
                    }
                };
                (Some(kind), mid, dir, len)
            }
            None => {
                let he = edge.half_edge;
                let p = brep
                    .half_edges
                    .get(he)
                    .and_then(|h| brep.vertices.get(h.origin))
                    .map(|v| v.position)
                    .unwrap_or(Point3::new(0.0, 0.0, 0.0));
                (None, p, None, 0.0)
            }
        };
        let he = edge.half_edge;
        let (faces, is_seam) = match brep.half_edges.get(he) {
            Some(h) => {
                let f1 = brep.loops.get(h.loop_).map(|l| l.face);
                let f2 = brep
                    .half_edges
                    .get(h.twin)
                    .and_then(|t| brep.loops.get(t.loop_))
                    .map(|l| l.face);
                let kind_of = |f: Option<popkern_core::topology::FaceHandle>| {
                    f.and_then(|fh| brep.faces.get(fh))
                        .and_then(|face| face.surface.as_ref())
                        .map(|s| SurfaceKind::of(&s.geometry))
                };
                let seam = matches!((f1, f2), (Some(a), Some(b)) if a == b);
                ((kind_of(f1), kind_of(f2)), seam)
            }
            None => ((None, None), false),
        };
        let endpoints = brep.half_edges.get(edge.half_edge).and_then(|h| {
            let a = brep.vertices.get(h.origin)?.position;
            let b = brep
                .half_edges
                .get(h.twin)
                .and_then(|t| brep.vertices.get(t.origin))?
                .position;
            Some((a, b))
        });
        out.push(EdgeInfo {
            endpoints,
            handle: eh,
            curve,
            midpoint,
            direction,
            faces,
            is_seam,
            chord_length,
            name: edge.name.and_then(|nid| brep.name_pool.get(nid).cloned()),
        });
    }
    out
}

/// Collect the queryable facts about every face of a body.
pub fn face_infos(brep: &BRep) -> Vec<FaceInfo> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let mut out = Vec::new();
    for (fh, face) in brep.faces.iter() {
        let kind = face.surface.as_ref().map(|s| SurfaceKind::of(&s.geometry));
        // Outward normal from the outer-loop WINDING (Newell), not
        // the chart normal ⊕ orientation flag: loop winding is the
        // representation-independent outwardness (measured: a sweep
        // cap's chart normal matched +Z on both caps).
        let mut acc = Vec3::new(0.0, 0.0, 0.0);
        let mut newell = Vec3::new(0.0, 0.0, 0.0);
        let mut count = 0usize;
        if let Some(l) = brep.loops.get(face.outer_loop) {
            let start = l.half_edge;
            let mut he = start;
            while let Some(d) = brep.half_edges.get(he) {
                let Some(next) = brep.half_edges.get(d.next) else {
                    break;
                };
                if let (Some(v), Some(vn)) =
                    (brep.vertices.get(d.origin), brep.vertices.get(next.origin))
                {
                    let (p, q) = (v.position, vn.position);
                    acc = acc + (p - Point3::new(0.0, 0.0, 0.0));
                    newell.x += (p.y - q.y) * (p.z + q.z);
                    newell.y += (p.z - q.z) * (p.x + q.x);
                    newell.z += (p.x - q.x) * (p.y + q.y);
                    count += 1;
                }
                he = d.next;
                if he == start {
                    break;
                }
            }
        }
        let normal = matches!(
            face.surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Plane { .. })
        )
        .then(|| newell.normalized())
        .flatten();
        let centroid = Point3::new(0.0, 0.0, 0.0) + acc * (1.0 / count.max(1) as f64);
        out.push(FaceInfo {
            handle: fh,
            kind,
            centroid,
            normal,
            name: face.name.and_then(|nid| brep.name_pool.get(nid).cloned()),
        });
    }
    out
}
