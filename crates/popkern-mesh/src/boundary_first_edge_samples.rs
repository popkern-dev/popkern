//! Canonical per-edge sampling for the boundary-first tessellator:
//! every edge is sampled ONCE, adaptively by 3D chord height, and
//! both incident faces consume the same polyline (watertight by
//! construction). Split out of `boundary_first.rs` (500-line cap).

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle};

/// Canonical per-edge samples: parameter fractions (0..=1 of the
/// edge's curve window) and their 3D points. Both faces adjacent to
/// the edge consume the same set.
pub(crate) struct EdgeSamples {
    /// Curve-window fractions, ascending, including 0 and 1.
    pub fractions: Vec<f64>,
    /// 3D points at those fractions.
    pub points: Vec<Point3>,
}

/// Sample every edge with a finite curve window once, adaptively by
/// 3D chord height.
pub(crate) fn sample_all_edges(brep: &BRep, chord_height: f64) -> HashMap<EdgeHandle, EdgeSamples> {
    let mut map = HashMap::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            // Curveless edge: synthesize vertex-chord samples so the
            // ring builder doesn't bail (a boolean can leave
            // degenerate zero-length curveless edges at discarded
            // tool corners; a missing entry made every face touching
            // one fall back to the cracking uniform-grid path). Zero-length
            // edges get NO entry — the ring builder skips them.
            let Some((va, vb)) = edge_endpoint_positions(brep, edge.half_edge) else {
                continue;
            };
            if va.distance_to(&vb) > 1.0e-9 {
                map.insert(
                    eh,
                    EdgeSamples {
                        fractions: vec![0.0, 1.0],
                        points: vec![va, vb],
                    },
                );
            }
            continue;
        };
        if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
            continue;
        }
        let eval: &dyn CurveEval = &curve.geometry;
        let at = |f: f64| eval.point_at(curve.t_start + (curve.t_end - curve.t_start) * f);
        // Adaptive subdivision by chord height, depth-limited.
        let mut fractions = vec![0.0, 1.0];
        let mut i = 0;
        while i + 1 < fractions.len() {
            let (a, b) = (fractions[i], fractions[i + 1]);
            let (pa, pm, pb) = (at(a), at(0.5 * (a + b)), at(b));
            let chord_mid = Point3::new(
                0.5 * (pa.x + pb.x),
                0.5 * (pa.y + pb.y),
                0.5 * (pa.z + pb.z),
            );
            if pm.distance_to(&chord_mid) > chord_height && fractions.len() < 256 && (b - a) > 1e-3
            {
                fractions.insert(i + 1, 0.5 * (a + b));
            } else {
                i += 1;
            }
        }
        let mut points: Vec<Point3> = fractions.iter().map(|&f| at(f)).collect();
        // Snap the endpoint samples to the shared VERTEX positions:
        // curve-endpoint evaluations can differ from the topological
        // vertex by ~1e-8 (refit/snap noise), and two edges meeting
        // at a vertex would then produce near-miss boundary points
        // that straddle the weld rounding. The gate is the body's own
        // linear tolerance (floored at 1e-6): a curve end within the
        // tolerance of its vertex IS that vertex by the body's
        // contract (a bored sphere at the 1e-3 tolerance carries its
        // exact seam arc's end 1e-4 from the vertex the 1e-3
        // intersection placed; an absolute 1e-6 gate left the seam
        // and rim rings disagreeing there and 16 open mesh edges).
        let snap_gate = brep.tolerance.linear().max(1.0e-6);
        if let (Some(first), Some(last), Some((va, vb))) = (
            points.first().copied(),
            points.last().copied(),
            edge_endpoint_positions(brep, edge.half_edge),
        ) {
            let (pa, pb) = if first.distance_to(&va) + last.distance_to(&vb)
                <= first.distance_to(&vb) + last.distance_to(&va)
            {
                (va, vb)
            } else {
                (vb, va)
            };
            if first.distance_to(&pa) < snap_gate {
                points[0] = pa;
            }
            let n = points.len() - 1;
            if points[n].distance_to(&pb) < snap_gate {
                points[n] = pb;
            }
        }
        map.insert(eh, EdgeSamples { fractions, points });
    }
    map
}

/// Endpoint vertex positions of an edge via its canonical half-edge.
///
/// Boundary edges have a dead twin; the end vertex then comes from the
/// surviving side's `next` origin. Returns `None` when neither
/// resolution succeeds (stale canonical half-edge) — callers skip the
/// edge rather than panic on partially-retired topology.
pub(crate) fn edge_endpoint_positions(
    brep: &BRep,
    half_edge: popkern_core::topology::HalfEdgeHandle,
) -> Option<(Point3, Point3)> {
    let he = brep.half_edges.get(half_edge)?;
    let end_vertex = match brep.half_edges.get(he.twin) {
        Some(twin) => twin.origin,
        None => brep.half_edges.get(he.next)?.origin,
    };
    Some((
        brep.vertices.get(he.origin)?.position,
        brep.vertices.get(end_vertex)?.position,
    ))
}

/// Env-gated Delaunay validity audit (the delaunay-audit trace topic): count
/// circumcircle violations in a Bowyer-Watson output. Proved the
/// PK-146 seam-band fans were LEGAL (zero violations) — the defect
/// was the cocircular lattice, not the triangulation.
pub(crate) fn delaunay_audit(
    face_handle: popkern_core::topology::FaceHandle,
    bw_pts: &[(f64, f64)],
    tris: &[[usize; 3]],
) {
    if !popkern_core::trace::trace_enabled("delaunay-audit") || bw_pts.len() <= 500 {
        return;
    }
    {
        let mut violations = 0usize;
        let mut worst = 0.0_f64;
        for t in tris {
            let (a, b, c) = (bw_pts[t[0]], bw_pts[t[1]], bw_pts[t[2]]);
            let d2 = 2.0 * (a.0 * (b.1 - c.1) + b.0 * (c.1 - a.1) + c.0 * (a.1 - b.1));
            if d2.abs() < 1e-30 {
                continue;
            }
            let ux = ((a.0 * a.0 + a.1 * a.1) * (b.1 - c.1)
                + (b.0 * b.0 + b.1 * b.1) * (c.1 - a.1)
                + (c.0 * c.0 + c.1 * c.1) * (a.1 - b.1))
                / d2;
            let uy = ((a.0 * a.0 + a.1 * a.1) * (c.0 - b.0)
                + (b.0 * b.0 + b.1 * b.1) * (a.0 - c.0)
                + (c.0 * c.0 + c.1 * c.1) * (b.0 - a.0))
                / d2;
            let r2 = (a.0 - ux).powi(2) + (a.1 - uy).powi(2);
            for (k, p) in bw_pts.iter().enumerate() {
                if k == t[0] || k == t[1] || k == t[2] {
                    continue;
                }
                let dd = (p.0 - ux).powi(2) + (p.1 - uy).powi(2);
                if dd < r2 * (1.0 - 1e-9) {
                    violations += 1;
                    let depth = (r2.sqrt() - dd.sqrt()) / r2.sqrt();
                    if depth > worst {
                        worst = depth;
                    }
                    break;
                }
            }
        }
        eprintln!(
            "  [bf-delaunay-audit] face {face_handle:?}: {} tris, {violations} circumcircle violations (worst depth {worst:.3})",
            tris.len()
        );
    }
}
