//! Periodic pcurve branch reconciliation.
//!
//! On a u-periodic surface every pcurve endpoint is only defined
//! modulo the period, so a loop can be endpoint-continuous while a
//! side's INTERNAL sweep runs the long way around the chart (its
//! pcurve Δu disagrees with the edge curve's true angular sweep by a
//! multiple of the period). Such a loop has net UV winding ±period
//! instead of 0, and every UV-domain consumer misreads the trim: the
//! measured face area integrates the full drum, tessellation covers
//! the complement, STEP re-import shifts the patch.
//!
//! `reconcile_periodic_pcurve_branches` repairs this per face:
//! 1. For each outer-loop side, the TRUE parameter sweep is measured
//!    branch-free by stepwise-continuation projection of edge-curve
//!    samples (consecutive samples are close, so nearest-branch
//!    unwrapping per step is unambiguous).
//! 2. A side whose pcurve sweep disagrees with the true sweep by more
//!    than half a period is rebuilt from the continuation samples
//!    (collapsed to an exact degree-1 line when the image is
//!    straight, a degree-1 polyline otherwise).
//! 3. Sides are then chained: each side's start branch is shifted by
//!    whole periods to continue the previous side's end.
//! 4. The result commits only if the reconciled loop closes with net
//!    zero u-winding; otherwise the face is left untouched.

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::topology::{FaceHandle, HalfEdgeHandle};

use super::brep::BRep;

impl BRep {
    /// Reconcile period branches of the outer-loop pcurves on a
    /// u-periodic face. No-op for non-periodic surfaces, faces
    /// without full pcurve coverage, and loops already at net zero
    /// u-winding. Returns true when the face's pcurves changed.
    pub fn reconcile_periodic_pcurve_branches(&mut self, face: FaceHandle) -> bool {
        let Some(face_data) = self.faces.get(face) else {
            return false;
        };
        let Some(surface) = face_data.surface.as_ref().map(|s| s.geometry.clone()) else {
            return false;
        };
        let Some(period) = surface.periodicity().0 else {
            return false;
        };
        let Some(loop_data) = self.loops.get(face_data.outer_loop) else {
            return false;
        };
        let mut sides = Vec::new();
        let start = loop_data.half_edge;
        let mut current = start;
        for _ in 0..self.half_edges.len() {
            sides.push(current);
            let Some(data) = self.half_edges.get(current) else {
                return false;
            };
            current = data.next;
            if current == start {
                break;
            }
        }
        if sides.len() < 2 {
            return false;
        }
        let mut spans = Vec::with_capacity(sides.len());
        for side in &sides {
            let Some(pcurve) = self.half_edges[*side].pcurve.as_ref() else {
                return false;
            };
            let (t0, t1) = pcurve.parameter_range();
            let s = pcurve.point_at(t0);
            let e = pcurve.point_at(t1);
            spans.push((s.x, e.x));
        }
        let net: f64 = spans.iter().map(|(s, e)| e - s).sum::<f64>()
            + spans
                .iter()
                .zip(spans.iter().cycle().skip(1))
                .map(|((_, e), (s_next, _))| nearest_branch_delta(s_next - e, period))
                .sum::<f64>();
        let trace = std::env::var_os("PK_PCURVE_TRACE").is_some();
        if trace {
            eprintln!(
                "[pcurve-period] face {face:?} net u-winding {net:.3} over {} side(s)",
                sides.len()
            );
        }
        if net.abs() < period / 4.0 {
            return false;
        }

        // Measure each side's true sweep and rebuild the outliers,
        // then chain every side's start branch to continue the
        // previous side's end. The mutation commits on a candidate
        // pcurve set and is applied only if the reconciled loop
        // closes at net zero u-winding.
        let mut candidate: Vec<Pcurve> = Vec::with_capacity(sides.len());
        for (index, side) in sides.iter().enumerate() {
            let Some(samples) = self.continuation_uv_samples(*side, &surface, period) else {
                if trace {
                    eprintln!("[pcurve-period]   side {side:?} continuation projection failed");
                }
                return false;
            };
            let true_sweep = samples.last().map(|p| p.x).unwrap_or(0.0)
                - samples.first().map(|p| p.x).unwrap_or(0.0);
            let (s_u, e_u) = spans[index];
            if trace {
                eprintln!(
                    "[pcurve-period]   side {side:?} pcurve sweep {:.3} true sweep {true_sweep:.3}",
                    e_u - s_u
                );
            }
            if (e_u - s_u - true_sweep).abs() > period / 2.0 {
                candidate.push(polyline_pcurve(&samples));
            } else {
                let Some(existing) = self.half_edges[*side].pcurve.clone() else {
                    return false;
                };
                candidate.push(existing);
            }
        }
        let mut previous_end: Option<f64> = None;
        for pcurve in candidate.iter_mut() {
            let (t0, t1) = pcurve.parameter_range();
            let s_u = pcurve.point_at(t0).x;
            if let Some(prev) = previous_end {
                let jump = s_u - prev;
                let shift = nearest_branch_delta(jump, period) - jump;
                if shift.abs() > period / 4.0 {
                    pcurve.shift_u(shift);
                }
            }
            previous_end = Some(pcurve.point_at(t1).x);
        }
        let reconciled_net: f64 = candidate
            .iter()
            .map(|pcurve| {
                let (t0, t1) = pcurve.parameter_range();
                pcurve.point_at(t1).x - pcurve.point_at(t0).x
            })
            .sum::<f64>()
            + candidate
                .iter()
                .zip(candidate.iter().cycle().skip(1))
                .map(|(a, b)| {
                    let (_, t1) = a.parameter_range();
                    let (t0, _) = b.parameter_range();
                    b.point_at(t0).x - a.point_at(t1).x
                })
                .sum::<f64>();
        if reconciled_net.abs() > period / 4.0 {
            if trace {
                eprintln!(
                    "[pcurve-period]   face {face:?} reconciled net {reconciled_net:.3} still \
                     wound — declined"
                );
            }
            return false;
        }
        for (side, pcurve) in sides.iter().zip(candidate) {
            self.half_edges[*side].pcurve = Some(pcurve);
        }
        true
    }

    /// Project edge-curve samples for one half-edge by stepwise
    /// continuation: each sample unwraps to the nearest branch of its
    /// predecessor, so the sweep is branch-free. Samples run in the
    /// half-edge's traversal direction.
    fn continuation_uv_samples(
        &self,
        side: HalfEdgeHandle,
        surface: &crate::geometry::surface::SurfaceGeometry,
        period: f64,
    ) -> Option<Vec<Point3>> {
        let data = self.half_edges.get(side)?;
        let curve = self.edges.get(data.edge)?.curve.as_ref()?;
        let origin = self.vertices.get(data.origin)?.position;
        let start_pos = curve.geometry.point_at(curve.t_start);
        let end_pos = curve.geometry.point_at(curve.t_end);
        let closed = start_pos.distance_to(&end_pos) < 1e-9;
        let forward = if closed {
            self.edges.get(data.edge)?.half_edge == side
        } else {
            origin.distance_to(&start_pos) <= origin.distance_to(&end_pos)
        };
        let (t_lo, t_hi) = if forward {
            (curve.t_start, curve.t_end)
        } else {
            (curve.t_end, curve.t_start)
        };
        const SAMPLES: usize = 17;
        let mut points = Vec::with_capacity(SAMPLES);
        let mut previous_u: Option<f64> = None;
        for i in 0..SAMPLES {
            let t = t_lo + (t_hi - t_lo) * i as f64 / (SAMPLES - 1) as f64;
            let p = curve.geometry.point_at(t);
            let (mut u, v) = project_point_to_uv(surface, &p)?;
            if let Some(prev) = previous_u {
                u = prev + nearest_branch_delta(u - prev, period);
            }
            previous_u = Some(u);
            points.push(Point3::new(u, v, 0.0));
        }
        Some(points)
    }
}

/// Shift `delta` by whole periods into (-period/2, period/2].
fn nearest_branch_delta(delta: f64, period: f64) -> f64 {
    let mut d = delta;
    while d > period / 2.0 {
        d -= period;
    }
    while d <= -period / 2.0 {
        d += period;
    }
    d
}

/// Degree-1 pcurve through the continuation samples; collapses to a
/// two-point exact line when the samples are straight in UV.
fn polyline_pcurve(samples: &[Point3]) -> Pcurve {
    let first = samples[0];
    let last = samples[samples.len() - 1];
    let chord = last - first;
    let len2 = chord.dot(&chord);
    let straight = samples.iter().all(|p| {
        if len2 < 1e-300 {
            return true;
        }
        let f = (*p - first).dot(&chord) / len2;
        let q = Point3::new(
            first.x + chord.x * f,
            first.y + chord.y * f,
            first.z + chord.z * f,
        );
        p.distance_to(&q) <= 1.0e-9
    });
    let control_points: Vec<Point3> = if straight {
        vec![first, last]
    } else {
        samples.to_vec()
    };
    let n = control_points.len();
    let mut knots = Vec::with_capacity(n + 2);
    knots.push(0.0);
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    Pcurve::from_nurbs(NurbsCurve {
        degree: 1,
        knots,
        control_points,
        weights: vec![1.0; n],
    })
}

#[cfg(test)]
#[path = "pcurve_period_branch_tests.rs"]
mod tests;
