//! Loop period winding: a face's trim must bound a region of its chart.
//!
//! On a chart periodic in u (a cylinder, a torus, a revolve) a loop can
//! wind around the period: a full cylinder band is bounded by one loop
//! winding +1 and one winding -1. The windings of all a face's loops
//! must cancel — the face is a region, so its boundary is null in the
//! chart's homology — unless the chart collapses a whole row to a
//! point (a sphere pole, a cone apex), where one winding loop bounds a
//! cap around the pole. A face whose loops leave a net winding has an
//! edge that runs the long way round the period (measured: a direct-
//! apex cap rim on a cube corner fillet ran three quarters round its
//! r = 0.1 strip, behind the axis through the solid's material; the
//! strip meshed the whole cylinder and the trim classifier's verdict on
//! it depended on where the loop polygon was sampled, PK-792). The
//! check reads each loop's pcurves only: within a pcurve the
//! displacement is continuous, and a join between consecutive pcurves
//! is taken to the nearest period image.

use crate::error::ValidationError;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, FaceHandle, LoopHandle};

/// How close the loop displacement must come to whole periods for the
/// winding to be read (a loop whose pcurve joins leave larger gaps is
/// someone else's defect: `OpenLoop`, check 9).
const WHOLE_PERIOD_SLACK: f64 = 0.25;

impl BRep {
    /// Net `(u, v)` period winding of `loop_handle` on `surface`, or
    /// `None` when a half-edge has no pcurve or the displacement is not
    /// near whole periods. An axis that is not periodic reads 0.
    pub fn loop_period_winding(
        &self,
        loop_handle: LoopHandle,
        surface: &SurfaceGeometry,
    ) -> Option<(i64, i64)> {
        let periods = surface.parameter_periodicity();
        let wrap = |d: f64, period: Option<f64>| match period {
            Some(p) if p > 0.0 => d - p * (d / p).round(),
            _ => d,
        };
        let start = self.loops.get(loop_handle)?.half_edge;
        let (mut du, mut dv) = (0.0, 0.0);
        let mut previous_end: Option<(f64, f64)> = None;
        let mut first_start: Option<(f64, f64)> = None;
        let mut current = start;
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(current)?;
            let pcurve = half_edge.pcurve.as_ref()?;
            let a = pcurve.point_at(pcurve.t_start);
            let b = pcurve.point_at(pcurve.t_end);
            if let Some(end) = previous_end {
                du += wrap(a.x - end.0, periods.0);
                dv += wrap(a.y - end.1, periods.1);
            }
            first_start.get_or_insert((a.x, a.y));
            du += b.x - a.x;
            dv += b.y - a.y;
            previous_end = Some((b.x, b.y));
            current = half_edge.next;
            if current == start {
                break;
            }
        }
        if current != start {
            return None;
        }
        let (first, last) = (first_start?, previous_end?);
        du += wrap(first.0 - last.0, periods.0);
        dv += wrap(first.1 - last.1, periods.1);
        let turns = |d: f64, period: Option<f64>| match period {
            Some(p) if p > 0.0 => {
                let k = d / p;
                ((k - k.round()).abs() <= WHOLE_PERIOD_SLACK).then_some(k.round() as i64)
            }
            _ => Some(0),
        };
        Some((turns(du, periods.0)?, turns(dv, periods.1)?))
    }

    /// Every face whose loops leave a net period winding on a chart
    /// with no collapsed row (see the module docs).
    pub(crate) fn validate_loop_period_winding_all(&self) -> Vec<ValidationError> {
        let faces: Vec<FaceHandle> = self.faces.keys().collect();
        self.validate_loop_period_winding_among(&faces)
    }

    /// [`ValidationError::LoopPeriodWinding`] for each of `faces` whose
    /// loops leave a net period winding on a chart with no collapsed
    /// row: an operator's region-scoped postcondition (a corner cap and
    /// the strips whose rims it rewrote).
    pub fn validate_loop_period_winding_among(&self, faces: &[FaceHandle]) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        for &face in faces {
            let Some(data) = self.faces.get(face) else {
                continue;
            };
            let Some(surface) = data.surface.as_ref().map(|s| &s.geometry) else {
                continue;
            };
            let (periodic_u, periodic_v) = surface.is_periodic();
            if !periodic_u && !periodic_v {
                continue;
            }
            let mut net = (0_i64, 0_i64);
            let mut determinate = true;
            for loop_handle in
                std::iter::once(data.outer_loop).chain(data.inner_loops.iter().copied())
            {
                match self.loop_period_winding(loop_handle, surface) {
                    Some((ku, kv)) => net = (net.0 + ku, net.1 + kv),
                    None => determinate = false,
                }
            }
            if !determinate || net == (0, 0) {
                continue;
            }
            let collapsed_u = net.0 != 0 && has_collapsed_row(surface, true);
            let collapsed_v = net.1 != 0 && has_collapsed_row(surface, false);
            if collapsed_u || collapsed_v {
                continue;
            }
            errors.push(ValidationError::LoopPeriodWinding {
                face,
                winding_u: net.0,
                winding_v: net.1,
            });
        }
        errors
    }
}

/// Whether some row across the periodic axis (`along_u`: a row of
/// constant v, varying u) collapses to one point, so a single winding
/// loop can bound a cap around it. Analytic charts answer by kind;
/// others test the rows at the finite ends of the other axis's range
/// (a revolve touching its axis), where a valid chart puts a collapse.
fn has_collapsed_row(surface: &SurfaceGeometry, along_u: bool) -> bool {
    match surface {
        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cone { .. } => return true,
        SurfaceGeometry::Cylinder { .. } | SurfaceGeometry::Torus { .. } => return false,
        _ => {}
    }
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let (fixed_range, vary_range) = if along_u {
        ((v0, v1), (u0, u1))
    } else {
        ((u0, u1), (v0, v1))
    };
    if !(vary_range.0.is_finite() && vary_range.1.is_finite()) {
        // No bounded row to test: make no claim.
        return true;
    }
    [fixed_range.0, fixed_range.1]
        .into_iter()
        .filter(|f| f.is_finite())
        .any(|fixed| {
            let at = |s: f64| {
                let t = vary_range.0 + (vary_range.1 - vary_range.0) * s;
                if along_u {
                    surface.point_at(t, fixed)
                } else {
                    surface.point_at(fixed, t)
                }
            };
            let (p0, p1, p2) = (at(0.0), at(1.0 / 3.0), at(2.0 / 3.0));
            let size = p0.distance_to(&at(0.5)).max(p1.distance_to(&p2));
            // A collapsed row: its thirds coincide while the chart
            // elsewhere has extent (scale from the opposite end row).
            let other = if (fixed - fixed_range.0).abs() < (fixed - fixed_range.1).abs() {
                fixed_range.1
            } else {
                fixed_range.0
            };
            let far = if other.is_finite() {
                let t = 0.5 * (vary_range.0 + vary_range.1);
                if along_u {
                    surface.point_at(t, other)
                } else {
                    surface.point_at(other, t)
                }
                .distance_to(&p0)
            } else {
                1.0
            };
            size <= 1.0e-9 * far.max(1.0e-12)
        })
}

#[cfg(test)]
#[path = "loop_period_winding_tests.rs"]
mod tests;
