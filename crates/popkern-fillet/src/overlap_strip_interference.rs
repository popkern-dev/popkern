//! PK-140 mid-face strip interference detection.
//!
//! Split from `overlap.rs` (file-size cap). See
//! [`detect_strip_interference`].

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;

/// A mid-face strip interference between two vertex-disjoint fillets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StripInterference {
    /// Spec index of the first fillet.
    pub fillet_a: usize,
    /// Spec index of the second fillet.
    pub fillet_b: usize,
}

/// Detect mid-face strip interference: two fillets on VERTEX-DISJOINT
/// edges sharing a support face whose removal strips overlap. Their
/// rails are typically PARALLEL (opposite edges of one face), so the
/// crossing test in [`detect_overlaps`] never fires — the overlap is
/// a region containment: one fillet's rail lies closer to the other
/// fillet's edge than that fillet's own rail does.
///
/// Corner pairs (edges sharing a vertex) are excluded — adjacent
/// fillets legitimately meet at vertex blends and their rails always
/// approach the shared corner.
pub fn detect_strip_interference(
    edge_geometries: &[crate::construct::EdgeGeometry],
    brep: &popkern_core::topology::BRep,
    tolerance: &Tolerance,
) -> Vec<StripInterference> {
    let samples = |curve: &NurbsCurve| -> Vec<popkern_core::geometry::Point3> {
        let (t0, t1) = curve.parameter_range();
        (0..=16)
            .map(|k| curve.point_at(t0 + (t1 - t0) * k as f64 / 16.0))
            .collect()
    };
    let edge_points = |index: usize| -> Vec<popkern_core::geometry::Point3> {
        let analysis = &edge_geometries[index].analysis;
        brep.edges
            .get(analysis.edge)
            .and_then(|edge| edge.curve.as_ref())
            .map(|curve| {
                (0..=16)
                    .map(|k| {
                        curve.geometry.point_at(
                            curve.t_start + (curve.t_end - curve.t_start) * k as f64 / 16.0,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default()
    };
    let min_distance =
        |a: &[popkern_core::geometry::Point3], b: &[popkern_core::geometry::Point3]| -> f64 {
            let mut best = f64::INFINITY;
            for pa in a {
                for pb in b {
                    best = best.min(pa.distance_to(pb));
                }
            }
            best
        };
    let mut interferences = Vec::new();
    for i in 0..edge_geometries.len() {
        for j in (i + 1)..edge_geometries.len() {
            let (a, b) = (&edge_geometries[i].analysis, &edge_geometries[j].analysis);
            // Corner pairs are legitimate vertex-blend territory.
            if a.vertex_start == b.vertex_start
                || a.vertex_start == b.vertex_end
                || a.vertex_end == b.vertex_start
                || a.vertex_end == b.vertex_end
            {
                continue;
            }
            // Shared support face pairs, with each fillet's trim on
            // that face.
            let shared: Vec<(&NurbsCurve, &NurbsCurve)> = [
                (a.face_1, &edge_geometries[i].trim_3d_1),
                (a.face_2, &edge_geometries[i].trim_3d_2),
            ]
            .into_iter()
            .flat_map(|(fa, trim_a)| {
                [
                    (b.face_1, &edge_geometries[j].trim_3d_1),
                    (b.face_2, &edge_geometries[j].trim_3d_2),
                ]
                .into_iter()
                .filter(move |(fb, _)| *fb == fa)
                .map(move |(_, trim_b)| (trim_a, trim_b))
            })
            .collect();
            if shared.is_empty() {
                continue;
            }
            let pts_edge_i = edge_points(i);
            let pts_edge_j = edge_points(j);
            if pts_edge_i.is_empty() || pts_edge_j.is_empty() {
                continue;
            }
            let hit = shared.iter().any(|(trim_a, trim_b)| {
                let rail_a = samples(trim_a);
                let rail_b = samples(trim_b);
                let width_a = min_distance(&rail_a, &pts_edge_i);
                let width_b = min_distance(&rail_b, &pts_edge_j);
                let b_into_a = min_distance(&rail_b, &pts_edge_i);
                let a_into_b = min_distance(&rail_a, &pts_edge_j);
                let slack = tolerance.linear().max(1.0e-9) * 10.0;
                b_into_a + slack < width_a || a_into_b + slack < width_b
            });
            if hit {
                interferences.push(StripInterference {
                    fillet_a: i,
                    fillet_b: j,
                });
            }
        }
    }
    interferences
}

/// PK-140 thin-wall precheck: returns `Some((requested_radius,
/// available_extent))` when a PLANAR support face's full extent
/// perpendicular to the filleted edge is smaller than the requested
/// radius — the rail would exit the face across its far boundary.
/// Curved supports are skipped (their reach is not a vertex-extent
/// question); the max-vertex-reach extent only errs permissive.
pub(crate) fn planar_support_extent_overrun(
    brep: &popkern_core::topology::BRep,
    analysis: &crate::analysis::EdgeAnalysis,
    radius: &crate::FilletRadius,
    tolerance: &Tolerance,
) -> Option<(f64, f64)> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let requested = match radius {
        crate::FilletRadius::Constant(r) => *r,
        crate::FilletRadius::Variable { start, end } => start.max(*end),
        crate::FilletRadius::ByPoints(points) => {
            points.iter().map(|(_, r)| *r).fold(0.0_f64, f64::max)
        }
    };
    if !requested.is_finite() || requested <= 0.0 {
        return None;
    }
    let edge_a = brep.vertices.get(analysis.vertex_start)?.position;
    let edge_b = brep.vertices.get(analysis.vertex_end)?.position;
    let edge_dir = (edge_b - edge_a).normalized()?;
    let slack = tolerance.linear().max(1.0e-9) * 10.0;
    for (face, surface) in [
        (analysis.face_1, &analysis.surface_1),
        (analysis.face_2, &analysis.surface_2),
    ] {
        if !matches!(surface, SurfaceGeometry::Plane { .. }) {
            continue;
        }
        let mut extent = 0.0_f64;
        for p in brep.face_outer_loop_positions(face) {
            let d = p - edge_a;
            let perp = d - edge_dir * d.dot(&edge_dir);
            extent = extent.max(perp.length());
        }
        if extent > 0.0 && requested > extent + slack {
            return Some((requested, extent));
        }
    }
    None
}
