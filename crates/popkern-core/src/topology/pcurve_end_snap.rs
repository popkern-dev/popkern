//! Whether endpoint surgery can repair a pcurve (PK-610).
//!
//! `fix_pcurve_endpoints` moves a NURBS pcurve's first or last control
//! point by the UV error at that end: sound for the drift a vertex weld
//! leaves, wrong when the pcurve spans a different stretch of the edge's
//! carrier. A setback carve shortened a cube strip's straight edge from
//! length 1.0 to 0.8; the strip's cubic pcurve still spanned the old
//! edge (its midpoint image lies on the shortened curve, so the
//! stale-image check kept it), and the surgery moved its last control
//! point by 0.2 in v against a last leg of 0.07. The result ended on
//! the right vertices but overshot the edge's end by 7.2e-2 in between:
//! a folded pcurve that validates (every image lies on the carrier
//! line) while trims and tessellation follow the fold.
//!
//! The verdict is the surgery's own postcondition: snap the ends as
//! `fix_pcurve_endpoints` would and require the result to image onto
//! the edge's trimmed curve. A drift repair passes; a mis-spanned
//! pcurve fails and is re-projected from the edge curve instead (a
//! half-the-end-leg shift bound was tried first: it caught the 0.2
//! shift but passed a smaller one that still bent the last span 3e-3
//! off a rim arc).

use super::HalfEdgeHandle;
use super::brep::BRep;
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Uniform samples (times four) when checking a snapped pcurve.
const SNAP_CHECK_SAMPLES: usize = 16;

impl BRep {
    /// Whether snapping `pcurve`'s ends to `heh`'s vertices (moving
    /// its first and last control points by the end UV errors) leaves
    /// it off the edge's trimmed curve by more than ten linear
    /// tolerances — the pcurve does not span this edge.
    pub(super) fn pcurve_end_snap_folds(
        &self,
        heh: HalfEdgeHandle,
        pcurve: &Pcurve,
        surface: &SurfaceGeometry,
    ) -> bool {
        use crate::geometry::curve::CurveEval;
        use crate::geometry::uv_projection::project_point_to_uv_near;
        let PcurveGeometry::Nurbs(nurbs) = &pcurve.geometry else {
            return false;
        };
        if nurbs.degree < 2 || nurbs.control_points.len() < 3 {
            return false;
        }
        let Some(he) = self.half_edges.get(heh) else {
            return false;
        };
        let Some(curve) = self.edges.get(he.edge).and_then(|e| e.curve.as_ref()) else {
            return false;
        };
        let Some(next) = self.half_edges.get(he.next) else {
            return false;
        };
        let (Some(origin), Some(dest)) =
            (self.vertices.get(he.origin), self.vertices.get(next.origin))
        else {
            return false;
        };
        let pos_tol = self.tolerance.linear().max(1e-9);
        // The snapped pcurve, exactly as endpoint surgery builds it.
        let mut snapped = nurbs.clone();
        let n = snapped.control_points.len();
        let mut moved = false;
        for (t, vertex, index) in [
            (pcurve.t_start, origin.position, 0),
            (pcurve.t_end, dest.position, n - 1),
        ] {
            let uv = pcurve.point_at(t);
            if surface.point_at(uv.x, uv.y).distance_to(&vertex) < pos_tol {
                continue;
            }
            let Some((u, v)) = project_point_to_uv_near(surface, &vertex, uv.x, uv.y) else {
                continue;
            };
            snapped.control_points[index].x += u - uv.x;
            snapped.control_points[index].y += v - uv.y;
            moved = true;
        }
        if !moved {
            return false;
        }
        let snapped = Pcurve {
            geometry: PcurveGeometry::Nurbs(snapped),
            t_start: pcurve.t_start,
            t_end: pcurve.t_end,
        };
        // Distance from a point to the TRIMMED edge curve: nearest
        // carrier parameter clamped to the edge's range.
        let geometry: &dyn CurveEval = &curve.geometry;
        let (lo, hi) = (
            curve.t_start.min(curve.t_end),
            curve.t_start.max(curve.t_end),
        );
        let off_curve = |p: &crate::geometry::Point3| -> f64 {
            geometry
                .closest_point(p)
                .map(|t| geometry.point_at(t.clamp(lo, hi)).distance_to(p))
                // Unmeasurable: keep the legacy surgery.
                .unwrap_or(0.0)
        };
        let gate = 10.0 * pos_tol;
        let span = pcurve.t_end - pcurve.t_start;
        // Uniform samples plus every knot span's midpoint: a snap moves
        // only the end spans, which may be a 1/128 sliver of a long
        // polyline pcurve.
        let uniform = (0..=4 * SNAP_CHECK_SAMPLES)
            .map(|k| pcurve.t_start + span * k as f64 / (4 * SNAP_CHECK_SAMPLES) as f64);
        let knots = &nurbs.knots;
        let midpoints = knots
            .windows(2)
            .filter(|w| w[1] > w[0])
            .map(|w| 0.5 * (w[0] + w[1]))
            .filter(|t| {
                (pcurve.t_start.min(pcurve.t_end)..=pcurve.t_start.max(pcurve.t_end)).contains(t)
            });
        uniform.chain(midpoints).any(|t| {
            let uv = snapped.point_at(t);
            off_curve(&surface.point_at(uv.x, uv.y)) > gate
        })
    }
}
