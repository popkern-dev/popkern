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
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
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

/// Whether chart points `a` and `b` lie on one pole: a line of the
/// chart along u that images one point (`S_u` vanishes there against
/// `S_v`), both points imaging it. Shifting a pcurve by a period to
/// close a u-gap on a pole moved a sphere face's second seam side onto
/// the first, and the face's loop then enclosed the strip between its
/// trims and the pole instead of the face (Viviani's cut with the ball's
/// pole at the double point meshed 3.0 short and measured 33.7 low,
/// PK-1139).
fn on_one_pole(surface: &SurfaceGeometry, a: (f64, f64), b: (f64, f64)) -> bool {
    let degenerate = |(u, v): (f64, f64)| {
        surface.partial_u(u, v).length() <= POLE_RATIO * surface.partial_v(u, v).length()
    };
    degenerate(a)
        && degenerate(b)
        && surface
            .point_at(a.0, a.1)
            .distance_to(&surface.point_at(b.0, b.1))
            <= POLE_RATIO * surface.partial_v(a.0, a.1).length()
}

/// The u-shift for `pc`, entering the loop across a pole: a whole
/// period when it is the second side of a seam whose first side the loop
/// already ran (its edge used earlier, its chart midpoint on that side's
/// branch), else none. A seam's two sides lie a period apart on the
/// chart; projection recompute cannot pick the branch, and where a
/// dimple split a ball's seam both sides came back on u = 0: the loop
/// then enclosed two slivers beside the seam instead of the face, and
/// the classify ray rejected every hit on the ball as off its trim (a
/// golf ball's last dimple kept its outer cap, PK-1194). The gap free
/// along the pole is spent where the loop needs it.
fn seam_side_across_pole(
    pc: &Pcurve,
    edge: super::EdgeHandle,
    first_use: &[(super::EdgeHandle, (f64, f64))],
    period: f64,
) -> f64 {
    let Some(&(_, (fu, fv))) = first_use.iter().find(|(e, _)| *e == edge) else {
        return 0.0;
    };
    let m = pc.point_at(0.5 * (pc.t_start + pc.t_end));
    let near = 1e-6 * period;
    if (m.x - fu).abs() > near || (m.y - fv).abs() > near {
        return 0.0;
    }
    // Toward the other border: up from the chart's lower half, down
    // from its upper (a NURBS chart re-bases the loop afterwards).
    if fu - period * (fu / period).floor() < 0.5 * period {
        period
    } else {
        -period
    }
}

/// `S_u` against `S_v` at a pole (a sphere's reads ~1e-16), and the
/// coincidence of the two points relative to the chart's scale.
const POLE_RATIO: f64 = 1e-9;

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
                // Each edge's first use in this loop: its chart midpoint
                // once shifted, for the seam-side test at a pole.
                let mut first_use: Vec<(super::EdgeHandle, (f64, f64))> = Vec::new();
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
                            // A gap along a pole line is no gap: every
                            // u there is one point, and the loop moves
                            // along it freely (a sphere face's seam
                            // sides meet the pole a period apart).
                            if !on_one_pole(&surface.geometry, pe, (s.x, s.y)) {
                                du = wrap(pe.0 - s.x, period_u);
                            } else if let Some(p) = period_u {
                                du = seam_side_across_pole(pc, h.edge, &first_use, p);
                            }
                            dv = wrap(pe.1 - s.y, period_v);
                        }
                        if du != 0.0 || dv != 0.0 {
                            shifts.push((he, du, dv));
                        }
                        if !first_use.iter().any(|(edge, _)| *edge == h.edge) {
                            let m = pc.point_at(0.5 * (pc.t_start + pc.t_end));
                            first_use.push((h.edge, (m.x + du, m.y + dv)));
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
                        .map(|s| s.geometry.parameter_range());
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
