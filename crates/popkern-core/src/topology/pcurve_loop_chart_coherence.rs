//! Per-loop chart coherence on periodic surface charts: branch
//! chaining, net-winding closure, and principal-cell normalization.
//!
//! The recompute path projects each half-edge independently, and cold
//! projections land in the principal branch — so two loop-mates can
//! come back a whole period apart even when both are individually
//! exact. The chart-valid invariant is per LOOP: after branch-chaining
//! each pcurve to its predecessor, the net winding must be zero and
//! the loop must sit in the surface's principal parameter cell.
//!
//! PK-286 made this invariant structural (enforced by
//! `recompute_face_pcurves` itself, not opt-in at batch wrappers) and
//! generalized it to BOTH periodic directions: the branch-honest
//! quarter-point chain aligns v on v-periodic charts (torus), so v
//! needs the same closure u always did.

use super::super::FaceHandle;
use super::super::brep::BRep;
use crate::geometry::pcurve::Pcurve;

/// Which chart coordinate a coherence pass operates on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChartAxis {
    /// The `u` coordinate (`Point3::x` of a pcurve evaluation).
    U,
    /// The `v` coordinate (`Point3::y` of a pcurve evaluation).
    V,
}

impl ChartAxis {
    fn coord(self, uv: crate::geometry::Point3) -> f64 {
        match self {
            ChartAxis::U => uv.x,
            ChartAxis::V => uv.y,
        }
    }

    fn shift(self, pcurve: &Pcurve, delta: f64) -> Pcurve {
        match self {
            ChartAxis::U => pcurve.translated_uv(delta, 0.0),
            ChartAxis::V => pcurve.translated_uv(0.0, delta),
        }
    }

    /// Magnitude of the surface partial along this axis at `uv` —
    /// zero at a parametric degeneracy of the axis (cone apex,
    /// sphere pole), where the axis coordinate carries no geometry.
    fn partial_len(
        self,
        surface: &crate::geometry::surface::SurfaceGeometry,
        uv: crate::geometry::Point3,
    ) -> f64 {
        use crate::geometry::surface::SurfaceEval;
        match self {
            ChartAxis::U => surface.partial_u(uv.x, uv.y).length(),
            ChartAxis::V => surface.partial_v(uv.x, uv.y).length(),
        }
    }
}

impl BRep {
    /// Zero the net winding of every loop of a face along every
    /// periodic chart direction, then normalize each loop to the
    /// principal parameter cell.
    ///
    /// The closed-edge pcurve direction convention (canonical
    /// half-edge forward, twin backward) makes the two uses of ONE
    /// edge consistent, but says nothing about a LOOP that chains
    /// several closed rims: a recess wall bounded by its bottom
    /// circle and the mouth ellipse got both rims wound the same
    /// way, the ring's net winding came to ±2 periods, and the
    /// boundary-first mesher rightly declined to chart it.
    /// Full-period members are flipped (via [`Pcurve::reversed`])
    /// until the net winding is zero.
    pub(in crate::topology) fn close_periodic_loop_winding(&mut self, face_handle: FaceHandle) {
        let Some((period_u, period_v)) = self
            .faces
            .get(face_handle)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.periodicity())
        else {
            return;
        };
        if let Some(period) = period_u.filter(|p| *p > 0.0) {
            self.close_loop_winding_along_axis(face_handle, period, ChartAxis::U);
        }
        if let Some(period) = period_v.filter(|p| *p > 0.0) {
            self.close_loop_winding_along_axis(face_handle, period, ChartAxis::V);
        }
    }

    fn close_loop_winding_along_axis(
        &mut self,
        face_handle: FaceHandle,
        period: f64,
        axis: ChartAxis,
    ) {
        let Some(face) = self.faces.get(face_handle) else {
            return;
        };
        let Some(surface_geom) = face.surface.as_ref().map(|s| s.geometry.clone()) else {
            return;
        };
        let mut loop_handles = vec![face.outer_loop];
        loop_handles.extend(face.inner_loops.iter().copied());

        for lh in loop_handles {
            let Some(loop_data) = self.loops.get(lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            if !self.half_edges.contains_key(start) {
                continue;
            }
            let mut hes = Vec::new();
            let mut cursor = start;
            let limit = self.half_edges.len() + 1;
            for _ in 0..limit {
                hes.push(cursor);
                cursor = self.half_edges[cursor].next;
                if cursor == start {
                    break;
                }
            }
            if hes.len() < 2 {
                continue;
            }

            // Iterate: chain branches, measure net winding, flip one
            // reversible full-period member, repeat. Bounded by the
            // loop length.
            for _ in 0..hes.len() {
                let mut prev_end: Option<f64> = None;
                let mut first_start: Option<f64> = None;
                let mut deltas: Vec<f64> = Vec::with_capacity(hes.len());
                let mut complete = true;
                for &heh in &hes {
                    let Some(pc) = self.half_edges[heh].pcurve.as_ref() else {
                        complete = false;
                        break;
                    };
                    let (t0, t1) = pc.parameter_range();
                    let joint_uv = pc.point_at(t0);
                    let mut sc = axis.coord(joint_uv);
                    let mut ec = axis.coord(pc.point_at(t1));
                    // A joint at a parametric degeneracy of this axis
                    // (cone apex, sphere pole: the axis partial
                    // vanishes) carries NO axis-continuity constraint —
                    // the legitimate branch jump of a periodic loop
                    // lives exactly there. Forcing k-shift continuity
                    // across it dragged an offset cone's seam pair onto
                    // one branch: the lateral loop tore at the seam,
                    // its trim polygon degenerated, and the mesh both
                    // cracked and lost the tail of the period (PK-332).
                    let mid_uv = pc.point_at(0.5 * (t0 + t1));
                    let joint_degenerate = axis.partial_len(&surface_geom, joint_uv)
                        < 1.0e-6 * (axis.partial_len(&surface_geom, mid_uv) + 1.0e-12);
                    if let Some(pe) = prev_end {
                        if joint_degenerate {
                            // Chain break: keep this member's authored
                            // branch and continue chaining from it.
                        } else {
                            let k = ((pe - sc) / period).round();
                            if k != 0.0 {
                                if let Some(pc_mut) = self.half_edges[heh].pcurve.as_mut() {
                                    *pc_mut = axis.shift(pc_mut, k * period);
                                }
                                sc += k * period;
                                ec += k * period;
                            }
                        }
                    } else {
                        first_start = Some(sc);
                    }
                    deltas.push(ec - sc);
                    prev_end = Some(ec);
                }
                if !complete {
                    break;
                }
                let (Some(fs), Some(pe)) = (first_start, prev_end) else {
                    break;
                };
                let net = pe - fs;
                if net.abs() <= 0.5 * period {
                    break;
                }
                let mut flipped = false;
                for (i, &heh) in hes.iter().enumerate() {
                    let d = deltas[i];
                    if d.abs() >= period * (1.0 - 1e-6)
                        && (net - 2.0 * d).abs() < net.abs() - 0.25 * period
                        && let Some(rev) =
                            self.half_edges[heh].pcurve.as_ref().map(|p| p.reversed())
                    {
                        self.half_edges[heh].pcurve = Some(rev);
                        flipped = true;
                        break;
                    }
                }
                if !flipped {
                    break;
                }
            }
            // PRINCIPAL-CELL NORMALIZATION: branch-chaining anchors
            // the loop to whatever cell its first pcurve landed on,
            // and analytic inversion of a rim whose arc starts past
            // +π puts the whole loop on the 2π..4π cell. That is
            // chart-equivalent for the kernel but OUTSIDE the
            // surface's principal parameter cell for interchange
            // (measured: OCC's BRepCheck rejected a twice-chamfered
            // column's bottom cone strip while the mesh welded
            // clean). Shift the WHOLE loop by whole periods so it
            // lands in [0, period) — a rigid translation that
            // preserves every relative relationship. (Every periodic
            // analytic chart anchors its periodic coordinate at 0,
            // so [0, period) IS the principal cell in both u and v.)
            {
                let mut min_c = f64::MAX;
                let mut max_c = f64::MIN;
                for &heh in &hes {
                    if let Some(pc) = self.half_edges[heh].pcurve.as_ref() {
                        let (t0, t1) = pc.parameter_range();
                        for t in [t0, t1] {
                            let c = axis.coord(pc.point_at(t));
                            min_c = min_c.min(c);
                            max_c = max_c.max(c);
                        }
                    }
                }
                // Scoped to loops WHOLLY outside [0, period): loops
                // that merely straddle a cell edge (charts anchored
                // at −π, seam-adjacent trims) are tuned behavior and
                // must not move (measured: shifting by the mean
                // broke eight direct-apex cap gates).
                if min_c <= max_c {
                    let eps = 1e-9 * period;
                    let k = if min_c >= period - eps {
                        (min_c / period + eps).floor()
                    } else if max_c <= eps {
                        (max_c / period - eps).ceil() - 1.0
                    } else {
                        0.0
                    };
                    if k != 0.0 {
                        for &heh in &hes {
                            if let Some(pc) = self.half_edges[heh].pcurve.as_mut() {
                                *pc = axis.shift(pc, -k * period);
                            }
                        }
                    }
                }
            }
        }
    }
}

impl BRep {
    /// Repair faces whose stored pcurve loops violate the periodic
    /// winding invariant by MORE than one period, by re-running
    /// [`BRep::recompute_face_pcurves`] (whose postcondition includes
    /// the per-loop closure) on exactly those faces. Returns the
    /// number of faces repaired.
    ///
    /// A net winding of ±1 period is the honest form for an annulus
    /// loop on a seamless periodic chart and is NOT touched; ≥2
    /// periods can only come from a loop whose pcurves were written
    /// by topology surgery AFTER the last recompute (measured: the
    /// flagship boss plate's mitred seam sliver mixed late-surgery
    /// seam pcurves with older rim pcurves, the boundary-first mesher
    /// declined the ±2-period ring, and the fallback grid left 98
    /// open weld edges). Pipelines that re-topologize faces after
    /// their last pcurve recompute call this as a final sweep.
    pub fn repair_periodic_loop_winding(&mut self) -> usize {
        let faces: Vec<FaceHandle> = self.faces.iter().map(|(fh, _)| fh).collect();
        let mut repaired = 0usize;
        for fh in faces {
            let Some((period_u, period_v)) = self
                .faces
                .get(fh)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.periodicity())
            else {
                continue;
            };
            let mut broken = false;
            for (axis, period) in [(ChartAxis::U, period_u), (ChartAxis::V, period_v)] {
                let Some(period) = period.filter(|p| *p > 0.0) else {
                    continue;
                };
                if self.loop_winding_exceeds(fh, period, axis, 1.5) {
                    broken = true;
                    break;
                }
            }
            if broken && self.recompute_face_pcurves(fh).is_ok() {
                repaired += 1;
            }
        }
        repaired
    }

    /// True when any loop of `face_handle` has stored-chain net
    /// winding along `axis` exceeding `factor` periods.
    fn loop_winding_exceeds(
        &self,
        face_handle: FaceHandle,
        period: f64,
        axis: ChartAxis,
        factor: f64,
    ) -> bool {
        let Some(face) = self.faces.get(face_handle) else {
            return false;
        };
        let mut loop_handles = vec![face.outer_loop];
        loop_handles.extend(face.inner_loops.iter().copied());
        for lh in loop_handles {
            let Some(loop_data) = self.loops.get(lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            if !self.half_edges.contains_key(start) {
                continue;
            }
            let mut prev_end: Option<f64> = None;
            let mut first_start: Option<f64> = None;
            let mut complete = true;
            let mut cursor = start;
            let limit = self.half_edges.len() + 1;
            for _ in 0..limit {
                let Some(pc) = self.half_edges[cursor].pcurve.as_ref() else {
                    complete = false;
                    break;
                };
                let (t0, t1) = pc.parameter_range();
                let mut sc = axis.coord(pc.point_at(t0));
                let mut ec = axis.coord(pc.point_at(t1));
                if let Some(pe) = prev_end {
                    let k = ((pe - sc) / period).round();
                    sc += k * period;
                    ec += k * period;
                } else {
                    first_start = Some(sc);
                }
                prev_end = Some(ec);
                cursor = self.half_edges[cursor].next;
                if cursor == start {
                    break;
                }
            }
            if !complete {
                continue;
            }
            if let (Some(fs), Some(pe)) = (first_start, prev_end)
                && (pe - fs).abs() > factor * period
            {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
#[path = "pcurve_loop_chart_coherence_tests.rs"]
mod pcurve_loop_chart_coherence_tests;
