//! Loop chart-continuity normalization for faces on geometrically
//! closed charts (PK-441).
//!
//! Validity model: a face's trim loop must be CONTINUOUS in its
//! chart — consecutive half-edge pcurves meet, with seam-edge sides
//! landing on opposite chart borders. Assembly pcurve adoption and
//! projection-based recompute cannot pick a seam BRANCH (u = 0 and
//! u = 1 image the same 3D curve on a closed chart), and a loop that
//! carries both seam sides on ONE border has chart discontinuities
//! that every chart consumer inherits: the doubly-ground coil's wall
//! meshed with hundreds of open edges, exactly the defect the sweep's
//! birth authorship fix removed for single-boolean output.
//!
//! The repair is exact: pcurves shift by WHOLE periods only, which is
//! image-preserving on a geometrically closed chart.

use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use crate::topology::{BRep, FaceHandle};

/// Translate a pcurve by `(du, dv)` in chart space.
fn translate_pcurve(pc: &mut Pcurve, du: f64, dv: f64) {
    match &mut pc.geometry {
        PcurveGeometry::Linear { start, end } => {
            start.0 += du;
            start.1 += dv;
            end.0 += du;
            end.1 += dv;
        }
        PcurveGeometry::Elliptical { center, .. } => {
            center.0 += du;
            center.1 += dv;
        }
        PcurveGeometry::Sinusoidal {
            u_offset, v_const, ..
        } => {
            *u_offset += du;
            *v_const += dv;
        }
        PcurveGeometry::Nurbs(curve) => {
            for cp in &mut curve.control_points {
                cp.x += du;
                cp.y += dv;
            }
        }
        PcurveGeometry::CylinderArcViaRationalQuad {
            theta_a,
            theta_b,
            v_const,
        } => {
            *theta_a += du;
            *theta_b += du;
            *v_const += dv;
        }
    }
}

impl BRep {
    /// Make every loop of `faces` chart-continuous on geometrically
    /// closed charts by shifting pcurves whole periods where
    /// consecutive pcurves jump by more than half a period. Faces on
    /// open charts, loops without pcurves, and jumps below the half-
    /// period gate are untouched. Anchored at each loop's first
    /// pcurve, so a loop already continuous is a no-op.
    pub fn normalize_loop_chart_continuity(&mut self, faces: &[FaceHandle]) {
        for &fh in faces {
            let Some(face) = self.faces.get(fh) else {
                continue;
            };
            let Some(surface) = face.surface.as_ref() else {
                continue;
            };
            let (period_u, period_v) =
                crate::geometry::chart_periods::geometric_chart_periods(&surface.geometry);
            if period_u.is_none() && period_v.is_none() {
                continue;
            }
            let mut loops = vec![face.outer_loop];
            loops.extend(face.inner_loops.iter().copied());
            for lh in loops {
                let Some(lp) = self.loops.get(lh) else {
                    continue;
                };
                let start_he = lp.half_edge;
                let mut he = start_he;
                let mut prev_end: Option<(f64, f64)> = None;
                // (handle, du, dv) shifts to apply after the walk.
                let mut shifts: Vec<(super::HalfEdgeHandle, f64, f64)> = Vec::new();
                for _ in 0..self.half_edges.len() + 1 {
                    let Some(h) = self.half_edges.get(he) else {
                        break;
                    };
                    let next = h.next;
                    if let Some(pc) = h.pcurve.as_ref() {
                        let s = pc.point_at(pc.t_start);
                        let e = pc.point_at(pc.t_end);
                        let mut du = 0.0;
                        let mut dv = 0.0;
                        if let (Some(pe), true) = (prev_end, true) {
                            let wrap = |gap: f64, period: Option<f64>| -> f64 {
                                match period {
                                    Some(p) if p.is_finite() && p > 0.0 && gap.abs() > 0.5 * p => {
                                        (gap / p).round() * p
                                    }
                                    _ => 0.0,
                                }
                            };
                            du = wrap(pe.0 - s.x, period_u);
                            dv = wrap(pe.1 - s.y, period_v);
                        }
                        if du != 0.0 || dv != 0.0 {
                            shifts.push((he, du, dv));
                        }
                        prev_end = Some((e.x + du, e.y + dv));
                    } else {
                        prev_end = None;
                    }
                    he = next;
                    if he == start_he {
                        break;
                    }
                }
                for (heh, du, dv) in shifts {
                    if let Some(h) = self.half_edges.get_mut(heh)
                        && let Some(pc) = h.pcurve.as_mut()
                    {
                        translate_pcurve(pc, du, dv);
                    }
                }
                // Re-base: a NURBS chart CLAMPS evaluation outside its
                // knot window, so the (now continuous) loop must live
                // inside it — anchoring on a wrong-branch first pcurve
                // can leave the whole loop a period off. Analytic
                // periodic charts evaluate period-invariantly and need
                // no re-base. The common shift is a whole period,
                // image-preserving.
                if matches!(
                    self.faces[fh].surface.as_ref().map(|s| &s.geometry),
                    Some(crate::geometry::surface::SurfaceGeometry::Nurbs(_))
                ) {
                    let range = self
                        .faces
                        .get(fh)
                        .and_then(|f| f.surface.as_ref())
                        .map(|s| {
                            use crate::geometry::surface::SurfaceEval;
                            s.geometry.parameter_range()
                        });
                    if let Some(((cu0, cu1), (cv0, cv1))) = range {
                        let (mut lo_u, mut hi_u) = (f64::MAX, f64::MIN);
                        let (mut lo_v, mut hi_v) = (f64::MAX, f64::MIN);
                        let start_he2 = self.loops[lh].half_edge;
                        let mut he2 = start_he2;
                        for _ in 0..self.half_edges.len() + 1 {
                            let Some(h) = self.half_edges.get(he2) else {
                                break;
                            };
                            if let Some(pc) = h.pcurve.as_ref() {
                                for t in [pc.t_start, 0.5 * (pc.t_start + pc.t_end), pc.t_end] {
                                    let q = pc.point_at(t);
                                    lo_u = lo_u.min(q.x);
                                    hi_u = hi_u.max(q.x);
                                    lo_v = lo_v.min(q.y);
                                    hi_v = hi_v.max(q.y);
                                }
                            }
                            he2 = h.next;
                            if he2 == start_he2 {
                                break;
                            }
                        }
                        let rebase = |lo: f64, hi: f64, c0: f64, c1: f64, period: Option<f64>| {
                            let Some(p) = period else { return 0.0 };
                            if !(p.is_finite() && p > 0.0) || lo > hi {
                                return 0.0;
                            }
                            // Whole periods moving the loop's span
                            // toward the chart window's span.
                            let mid = 0.5 * (lo + hi);
                            let cmid = 0.5 * (c0 + c1);
                            ((cmid - mid) / p).round() * p
                        };
                        let du = rebase(lo_u, hi_u, cu0, cu1, period_u);
                        let dv = rebase(lo_v, hi_v, cv0, cv1, period_v);
                        if du != 0.0 || dv != 0.0 {
                            let mut he2 = start_he2;
                            for _ in 0..self.half_edges.len() + 1 {
                                let Some(h) = self.half_edges.get_mut(he2) else {
                                    break;
                                };
                                let next = h.next;
                                if let Some(pc) = h.pcurve.as_mut() {
                                    translate_pcurve(pc, du, dv);
                                }
                                he2 = next;
                                if he2 == start_he2 {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
#[path = "loop_chart_continuity_tests.rs"]
mod loop_chart_continuity_tests;
