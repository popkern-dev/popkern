//! Split-stage cell crossings for doubly-periodic faces (extracted
//! from `chart_cell.rs` under the file-size ratchet): boundary-run
//! attribution (`nearest_boundary`) and the dense-sample + bisection
//! crossing scan (`cell_crossings`).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::{ChartCell, cell_windows, face_chart_cell};

impl ChartCell {
    /// Locate the boundary half-edge nearest a folded UV point, with
    /// the boundary parameter at the closest approach. Used to
    /// attribute a window end (an in/out transition of the query) to
    /// the boundary edge it crosses, for the split stage's
    /// vertex-creation machinery.
    pub(crate) fn nearest_boundary(&self, fu: f64, fv: f64) -> Option<(HalfEdgeHandle, f64, f64)> {
        let mut best: Option<(HalfEdgeHandle, f64, f64)> = None;
        for run in &self.ring {
            let n = run.samples.len();
            for i in 0..n.saturating_sub(1) {
                let a = run.samples[i];
                let b = run.samples[i + 1];
                let (dx, dy) = (b.0 - a.0, b.1 - a.1);
                let l2 = dx * dx + dy * dy;
                let s = if l2 < 1e-30 {
                    0.0
                } else {
                    (((fu - a.0) * dx + (fv - a.1) * dy) / l2).clamp(0.0, 1.0)
                };
                let (px, py) = (a.0 + s * dx, a.1 + s * dy);
                let d2 = (fu - px) * (fu - px) + (fv - py) * (fv - py);
                if best.as_ref().is_none_or(|(_, _, bd)| d2 < *bd) {
                    let frac = (i as f64 + s) / (n - 1) as f64;
                    let t_boundary = run.t0 + (run.t1 - run.t0) * frac;
                    best = Some((run.he, t_boundary, d2));
                }
            }
        }
        best.map(|(he, t, d2)| (he, t, d2.sqrt()))
    }
}

/// Cell-based boundary crossings for the SPLIT stage: window ends of
/// the folded query, each attributed to the boundary half-edge it
/// crosses (PK-304 — the legacy crossing search silently missed
/// splits on negative-branch and seam-band cells, leaving disc
/// regions connected and classification flood-fill leaking through).
/// `None` when the face is not cell-resolvable.
pub(crate) fn cell_crossings(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &popkern_core::geometry::nurbs_curve::NurbsCurve,
    tolerance: &Tolerance,
) -> Option<Vec<crate::clip::clip_uv::UvCrossing>> {
    use popkern_core::geometry::curve::CurveEval as _;
    let on = std::env::var_os("POPKERN_CELL_CLIP")
        .map(|v| v != "0")
        .unwrap_or(true);
    if !on {
        return None;
    }
    let cell = face_chart_cell(brep, face, tolerance)?;
    let (t0, t1) = curve_uv.parameter_range();
    let windows = cell_windows(&cell, curve_uv);
    let mut crossings = Vec::new();
    let _ = (t0, t1);
    for (enter, exit) in windows {
        for t_cross in [enter, exit] {
            // Window ends at the CURVE ends are kept when they sit on
            // the boundary: pair-level clipping already trims records
            // to in-face spans, so a clipped record's endpoints ARE
            // boundary crossings — the split stage needs them to
            // create junction vertices (skipping them re-broke every
            // split: 0 crossings, midpoint-inside, no partition).
            // Interior truncations are rejected by the
            // boundary-distance gate below either way.
            let q = curve_uv.point_at(t_cross);
            let (fu, fv) = cell.fold(q.x, q.y);
            let (he, t_boundary, dist) = cell.nearest_boundary(fu, fv)?;
            // The transition must actually sit on the boundary — a
            // fold artifact far from the ring is not a crossing.
            let scale = (cell.u_range.1 - cell.u_range.0)
                .abs()
                .max((cell.v_range.1 - cell.v_range.0).abs());
            if dist > scale * 1.0e-2 {
                continue;
            }
            crossings.push(crate::clip::clip_uv::UvCrossing {
                t_curve: t_cross,
                boundary_he: he,
                t_boundary,
            });
        }
    }
    Some(crossings)
}
