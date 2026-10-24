//! Boundary-riding stretches in the chart-cell clip (PK-942).
//!
//! A record can run ALONG a cell's boundary: a box face through the
//! torus axis cuts the primitive torus exactly on its meridian patch
//! edges. Strict even-odd containment decides a point on the ring by
//! floating-point noise. Unrotated, the record and the edge read the
//! same exact u and the test is consistently half-open (a point on the
//! low-u edge is inside, on the high-u edge outside), so exactly one of
//! the edge's two faces claims the record. Under a rotation about the
//! axis the edge reads 0.7 and the record 0.7 ± 3e-14: the verdict
//! flipped sample by sample, both faces claimed ~40 shreds each, and
//! the merge stranded the record sides (measured: 11 of 16 rotation_z
//! turns declined and one served a body whose torus area was 5.3x off).
//!
//! A maximal run of two or more samples within the coincidence band of
//! the ring is coincidence, not crossing. The whole run takes the
//! half-open verdict of the edge it rides, judged
//! at a point nudged off the ring along that edge's interior normal:
//! claimed when the interior lies toward +u (toward +v on an edge
//! along u). The two faces sharing the edge see opposite normals, so
//! exactly one claims it, in every orientation; unrotated the verdict
//! equals the strict one. The verdict is geometric even where the
//! strict test happens to agree along the run: each face pair marches
//! its own copy of the record, and the copies' noise differs (measured,
//! mixed rotation: face 2's copy of one cut flipped and face 6's sat
//! consistently on face 2's side, so a flip-gated override dropped the
//! cut from both). Samples off the ring keep the strict test.
//!
//! **Tolerance:** the band is the model's linear tolerance read in
//! uv; a record within it of the boundary over consecutive samples
//! coincides with the boundary by the model's own standard.

use super::ChartCell;
use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

impl ChartCell {
    /// Distance from a folded point to the nearest ring segment, and
    /// that segment's unit interior normal.
    fn nearest_segment(&self, fu: f64, fv: f64) -> Option<(f64, (f64, f64))> {
        let mut best: Option<(f64, (f64, f64))> = None;
        for run in &self.ring {
            for w in run.samples.windows(2) {
                let (a, b) = (w[0], w[1]);
                let (dx, dy) = (b.0 - a.0, b.1 - a.1);
                let l2 = dx * dx + dy * dy;
                if l2 < 1e-30 {
                    continue;
                }
                let s = (((fu - a.0) * dx + (fv - a.1) * dy) / l2).clamp(0.0, 1.0);
                let d = ((fu - a.0 - s * dx).powi(2) + (fv - a.1 - s * dy).powi(2)).sqrt();
                if best.is_none_or(|(bd, _)| d < bd) {
                    let l = l2.sqrt();
                    // Interior lies left of a counter-clockwise ring.
                    let n = if self.counter_clockwise {
                        (-dy / l, dx / l)
                    } else {
                        (dy / l, -dx / l)
                    };
                    best = Some((d, n));
                }
            }
        }
        best
    }
}

/// Per-sample containment nudges: `Some(offset)` on the samples of a
/// riding run (see the module docs), the uv offset at which the run's
/// half-open verdict is read; `None` elsewhere.
pub(super) fn riding_nudges(
    cell: &ChartCell,
    curve_uv: &NurbsCurve,
    ts: &[f64],
) -> Vec<Option<(f64, f64)>> {
    let samples: Vec<Option<(f64, f64)>> = ts
        .iter()
        .map(|&t| {
            let q = curve_uv.point_at(t);
            let (fu, fv) = cell.fold(q.x, q.y);
            cell.nearest_segment(fu, fv)
                .and_then(|(d, n)| (d <= cell.band).then_some(n))
        })
        .collect();
    let mut nudges = vec![None; ts.len()];
    let mut i = 0;
    while i < samples.len() {
        if samples[i].is_none() {
            i += 1;
            continue;
        }
        let start = i;
        while i < samples.len() && samples[i].is_some() {
            i += 1;
        }
        let run = &samples[start..i];
        if run.len() < 2 {
            continue;
        }
        let Some((nu, nv)) = run[run.len() / 2] else {
            continue;
        };
        let claimed = if nu.abs() > 1e-9 { nu > 0.0 } else { nv > 0.0 };
        let step = if claimed { cell.band } else { -cell.band };
        for nudge in &mut nudges[start..i] {
            *nudge = Some((nu * step, nv * step));
        }
    }
    nudges
}

#[cfg(test)]
#[path = "chart_cell_riding_tests.rs"]
mod tests;
