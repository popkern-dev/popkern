//! Corner-cap boundary arc reconstruction.
//!
//! The rolling-ball corner cap's true boundary curves are the
//! adjacent fillet strips' END CROSS-SECTIONS: circular arcs of the
//! fillet radius centered at the cap sphere's own center. They are
//! TANGENCY curves — no surface-intersection dispatcher can produce
//! them (the cap is tangent to its strips along them, and tangent
//! configurations are degenerate for intersection) — so they must be
//! reconstructed from the cap sphere itself. The carve phase leaves
//! straight chords in their place (2.5e-1 off-surface at demo
//! scales), which poisons tessellation and export association for
//! every face at the corner.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;

/// For each corner patch whose face is a rolling-ball Sphere,
/// replace chord boundary edges whose endpoints lie (within the
/// anchor-placement class, 5e-3) on the sphere with exact circular
/// arcs centered at the sphere center. Endpoint vertices are snapped
/// radially onto the sphere and adjacent chord-line edges rebuilt.
/// Returns the number of edges rewritten.
pub(super) fn arcify_corner_cap_boundaries(
    brep: &mut BRep,
    corner_patches: &[crate::types::CornerPatchInfo],
) -> usize {
    let mut rewritten = 0;
    for cp in corner_patches {
        let Some(face) = brep.faces.get(cp.face) else {
            continue;
        };
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Sphere { center, radius } = sof.geometry else {
            continue;
        };
        let start_he = brep.loops[face.outer_loop].half_edge;
        // Collect the loop's edges first (mutation below).
        let mut edges = Vec::new();
        let mut he = start_he;
        for _ in 0..brep.half_edges.len() + 1 {
            let d = &brep.half_edges[he];
            edges.push(d.edge);
            he = d.next;
            if he == start_he {
                break;
            }
        }
        for eh in edges {
            let Some(edge) = brep.edges.get(eh) else {
                continue;
            };
            let Some(curve) = edge.curve.as_ref() else {
                continue;
            };
            // Only chords: straight Lines or 2-point degree-1 NURBS.
            let is_chord = match &curve.geometry {
                CurveGeometry::Line { .. } => true,
                CurveGeometry::Nurbs(n) => n.degree == 1 && n.control_points.len() == 2,
                _ => false,
            };
            if !is_chord {
                continue;
            }
            let he0 = edge.half_edge;
            let (v1, v2) = {
                let d = &brep.half_edges[he0];
                (d.origin, brep.half_edges[d.twin].origin)
            };
            if v1 == v2 {
                continue;
            }
            let p1 = brep.vertices[v1].position;
            let p2 = brep.vertices[v2].position;
            let (d1, d2) = (p1.distance_to(&center), p2.distance_to(&center));
            let snap_gate = 5.0e-3;
            if (d1 - radius).abs() > snap_gate || (d2 - radius).abs() > snap_gate {
                continue;
            }
            if d1 < 1e-12 || d2 < 1e-12 {
                continue;
            }
            // Snap endpoints radially onto the sphere.
            let s1 = center + (p1 - center) * (radius / d1);
            let s2 = center + (p2 - center) * (radius / d2);
            crate::construct::edge_reconcile::snap_vertex_and_adjacent_chords(brep, v1, s1, eh);
            crate::construct::edge_reconcile::snap_vertex_and_adjacent_chords(brep, v2, s2, eh);
            // Exact minor arc centered at the sphere center.
            let a = s1 - center;
            let b = s2 - center;
            let Some(normal) = a.cross(&b).normalized() else {
                continue;
            };
            let Some(e1) = a.normalized() else { continue };
            let e2 = normal.cross(&e1);
            let sweep = (b.dot(&e2)).atan2(b.dot(&e1));
            // atan2 ∈ (−π, π]; a positive sweep traces the minor arc
            // from s1 to s2 around `normal`.
            let (arc_normal, sweep) = if sweep >= 0.0 {
                (normal, sweep)
            } else {
                (normal * -1.0, -sweep)
            };
            // Start angle is 0 in the arc's own (e1, n×e1) frame; the
            // Arc variant measures angles in the frame derived from
            // its normal, so re-derive theta0 for s1 in that frame.
            let (f1, f2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(arc_normal);
            let theta1 = a.dot(&f2).atan2(a.dot(&f1));
            let theta2_raw = b.dot(&f2).atan2(b.dot(&f1));
            let mut theta2 = theta2_raw;
            while theta2 < theta1 {
                theta2 += std::f64::consts::TAU;
            }
            // Consistency: the sweep should match; if not (wrapped
            // the long way), bring it into (theta1, theta1 + 2π).
            if (theta2 - theta1) - sweep > 1e-6 && theta2 - theta1 > std::f64::consts::PI {
                // Long way — the minor arc runs the other direction;
                // flip the normal-frame.
                continue;
            }
            let new_curve = CurveOnEdge {
                geometry: CurveGeometry::Arc {
                    center,
                    normal: arc_normal,
                    radius,
                    start_angle: theta1,
                    end_angle: theta2,
                },
                t_start: theta1,
                t_end: theta2,
            };
            // Verify the arc's midpoint is genuinely nearer both
            // endpoints' faces than the chord was (sanity: minor
            // arc, on-sphere).
            let ceval: &dyn CurveEval = &new_curve.geometry;
            let mid = ceval.point_at(0.5 * (theta1 + theta2));
            if (mid.distance_to(&center) - radius).abs() > 1e-9 {
                continue;
            }
            if let Some(e) = brep.edges.get_mut(eh) {
                e.curve = Some(new_curve);
                rewritten += 1;
            }
        }
    }
    rewritten
}
