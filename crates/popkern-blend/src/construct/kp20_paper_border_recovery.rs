//! Owning-surface recovery for stale face-edge handles in the KP20
//! paper-border extractor (submodule of `kp20_paper_border_extract`).

use popkern_core::geometry::{NurbsCurve, Point3, SurfaceGeometry};
use popkern_core::topology::BRep;

use crate::blend::{BlendBoundary, BoundarySegment};

/// Recover the owning surface of a face-edge segment whose recorded
/// face handle no longer exists in the BRep (the n-gon boundary is
/// computed before the corner-face topology surgery). The owner is
/// identified by incidence: its surface must contain this segment's
/// interior samples, but must NOT span the rest of the boundary loop
/// — the blend-face surface interpolates every boundary segment, so
/// any candidate that also contains most of the OTHER face-edge
/// midpoints is the cap face, not a support face. The incidence
/// tolerance is relative to this segment's chord (a local,
/// per-segment scale): on-face samples land at projection accuracy
/// (<< 1e-6 of chord), while a wrong face is off by a geometric
/// fraction of the chord.
pub(super) fn recover_face_edge_surface(
    brep: &BRep,
    curve: &NurbsCurve,
    boundary: &BlendBoundary,
) -> Result<SurfaceGeometry, String> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let point_on = |c: &NurbsCurve, t: f64| -> Point3 {
        let (t0, t1) = c.parameter_range();
        c.point_at(t0 + t * (t1 - t0))
    };
    let samples = [
        point_on(curve, 0.25),
        point_on(curve, 0.5),
        point_on(curve, 0.75),
    ];
    let chord = point_on(curve, 0.0).distance_to(&point_on(curve, 1.0));
    let tolerance = 1.0e-6 * chord.max(1.0e-9);
    let deviation = |surface: &SurfaceGeometry, p: &Point3| -> f64 {
        match project_point_to_uv(surface, p) {
            Some((u, v)) => surface.point_at(u, v).distance_to(p),
            None => f64::INFINITY,
        }
    };
    // Midpoints of the OTHER face-edge segments: a support face owns
    // one boundary segment, the blend-face surface owns them all.
    let foreign_midpoints: Vec<Point3> = boundary
        .segments
        .iter()
        .filter_map(|segment| match segment {
            BoundarySegment::FaceEdge { curve: other, .. } if !std::ptr::eq(other, curve) => {
                Some(point_on(other, 0.5))
            }
            _ => None,
        })
        .collect();
    let mut best: Option<(f64, SurfaceGeometry)> = None;
    for (_, face) in brep.faces.iter() {
        let Some(surface) = face.surface.as_ref().map(|s| &s.geometry) else {
            continue;
        };
        let worst = samples
            .iter()
            .map(|p| deviation(surface, p))
            .fold(0.0_f64, f64::max);
        if worst > tolerance {
            continue;
        }
        let foreign_hits = foreign_midpoints
            .iter()
            .filter(|p| deviation(surface, p) <= tolerance)
            .count();
        if 2 * foreign_hits > foreign_midpoints.len() {
            continue; // spans the loop: the blend-face surface, not the owner
        }
        if best.as_ref().is_none_or(|(d, _)| worst < *d) {
            best = Some((worst, surface.clone()));
        }
    }
    best.map(|(_, surface)| surface).ok_or_else(|| {
        format!(
            "stale face-edge handle: no live face surface contains the segment \
             (chord={chord:.3e}, tolerance={tolerance:.3e})"
        )
    })
}
