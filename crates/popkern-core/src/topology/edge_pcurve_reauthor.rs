//! Re-author an edge's two half-edge pcurves for a new spelling of its
//! curve, co-parameterised with that curve and running each half-edge
//! the way its existing pcurve runs.
//!
//! Used when an operator replaces an edge's geometry by an equivalent
//! form — the boolean's exact-rim promotion (NURBS circle → `Arc`,
//! PK-417) and the G2 fillet's containment the other way (`Arc` →
//! rational circle, PK-419). The loop's established UV sense is kept,
//! not the closed-rim recompute convention: a boolean rim's two faces
//! may both trace it clockwise (measured on a bored block), and a
//! convention-aligned re-author tore the seam-bounded wall's loop.
//! Pcurves are parameterised by the curve's own range; a half-edge
//! whose old trace runs against the curve gets the reversed pcurve
//! `P'(t) = P(t_start + t_end - t)`; periodic charts keep the old
//! pcurve's sheet; every pcurve is verified on the surface before it
//! is returned.

use crate::geometry::chart_pcurve::{analytic_chart_pcurve, fit_chart_pcurve};
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{NurbsCurve, Point3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

const VERIFY_SAMPLES: usize = 9;

/// Both half-edges' pcurves for `curve` as the new geometry of `edge`,
/// or `None` when a pcurve cannot be authored or verified on either
/// face (the caller keeps the old spelling). `gate` is the distance
/// each verification sample may miss the curve by.
pub fn coparameterised_pcurves(
    brep: &BRep,
    edge: EdgeHandle,
    curve: &CurveOnEdge,
    gate: f64,
    tolerance: &Tolerance,
) -> Option<Vec<(HalfEdgeHandle, Pcurve)>> {
    let e = brep.edges.get(edge)?;
    let forward = e.half_edge;
    let backward = brep.half_edges.get(forward)?.twin;
    let range = (curve.t_start, curve.t_end);
    let mut out = Vec::with_capacity(2);
    for (he, default_reversed) in [(forward, false), (backward, true)] {
        let half_edge = brep.half_edges.get(he)?;
        let face = brep.loops.get(half_edge.loop_)?.face;
        let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
        let ahead = pcurve_for(surface, curve, range, tolerance)?;
        let behind = reversed_pcurve(&ahead);
        let reversed = match half_edge.pcurve.as_ref() {
            Some(old) => {
                let trace = |p: &Pcurve| -> Vec<Point3> {
                    [0.25, 0.5, 0.75]
                        .iter()
                        .map(|f| {
                            let uv = p.point_at(p.t_start + (p.t_end - p.t_start) * f);
                            surface.point_at(uv.x, uv.y)
                        })
                        .collect()
                };
                let old_trace = trace(old);
                let score = |p: &Pcurve| -> f64 {
                    trace(p)
                        .iter()
                        .zip(&old_trace)
                        .map(|(a, b)| a.distance_to(b))
                        .fold(0.0, f64::max)
                };
                score(&behind) < score(&ahead)
            }
            None => default_reversed,
        };
        let mut pcurve = if reversed { behind } else { ahead };
        if let Some(old) = half_edge.pcurve.as_ref() {
            align_period(&mut pcurve, old, surface);
        }
        if !verified(surface, &pcurve, &curve.geometry, range, reversed, gate) {
            return None;
        }
        out.push((he, pcurve));
    }
    Some(out)
}

/// The pcurve of `curve` on `surface`: the exact plane arm, the
/// iso-line of a chart circle on a cylinder or cone, else the
/// co-parameterised chart fit.
fn pcurve_for(
    surface: &SurfaceGeometry,
    curve: &CurveOnEdge,
    range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Pcurve> {
    if let Some(exact) = analytic_chart_pcurve(surface, &curve.geometry, range) {
        return Some(exact);
    }
    if let Some(iso) = iso_circle_pcurve(surface, &curve.geometry, range) {
        return Some(iso);
    }
    let fitted = fit_chart_pcurve(surface, &curve.geometry, range, tolerance)?;
    Some(Pcurve {
        geometry: PcurveGeometry::Nurbs(fitted),
        t_start: range.0,
        t_end: range.1,
    })
}

/// A rim that is a chart iso-circle of a cylinder or cone (its plane
/// perpendicular to the axis, its centre on the axis) gets the
/// iso-linear pcurve the primitive and sweep builders spell — `u` the
/// chart angle, `v` the station — parameterised by the arc's angles.
fn iso_circle_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveGeometry,
    range: (f64, f64),
) -> Option<Pcurve> {
    let CurveGeometry::Arc {
        center,
        normal,
        radius,
        ..
    } = curve
    else {
        return None;
    };
    let (origin, axis, station, chart_radius) = match surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius: r,
        } => {
            let v = (*center - *origin).dot(axis);
            (*origin, *axis, v, *r)
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let v = (*center - *apex).dot(axis);
            (*apex, *axis, v, (v * half_angle.tan()).abs())
        }
        _ => return None,
    };
    let axial = axis.normalized()?;
    let aligned = normal.dot(&axial).abs() > 1.0 - 1.0e-9;
    let on_axis = {
        let d = *center - origin;
        (d - axial * d.dot(&axial)).length() < 1.0e-9 * radius.max(1.0)
    };
    if !(aligned && on_axis && (chart_radius - radius).abs() < 1.0e-9 * radius.max(1.0)) {
        return None;
    }
    let (e1, e2) = crate::geometry::analytic_surfaces::orthonormal_basis(axial);
    let start = curve.point_at(range.0) - origin;
    let u0 = start.dot(&e2).atan2(start.dot(&e1));
    let sense = if normal.dot(&axial) >= 0.0 { 1.0 } else { -1.0 };
    let sweep = range.1 - range.0;
    Some(Pcurve {
        geometry: PcurveGeometry::Linear {
            start: (u0, station),
            end: (u0 + sense * sweep, station),
        },
        t_start: range.0,
        t_end: range.1,
    })
}

/// Shift `pcurve` by whole periods so its start lands on the sheet of
/// `old`'s start.
fn align_period(pcurve: &mut Pcurve, old: &Pcurve, surface: &SurfaceGeometry) {
    let (period_u, period_v) = surface.parameter_periodicity();
    let new_start = pcurve.point_at(pcurve.t_start);
    let old_start = old.point_at(old.t_start);
    let shift = |period: Option<f64>, new: f64, old: f64| -> f64 {
        match period {
            Some(p) if p > 0.0 => p * ((old - new) / p).round(),
            _ => 0.0,
        }
    };
    let (du, dv) = (
        shift(period_u, new_start.x, old_start.x),
        shift(period_v, new_start.y, old_start.y),
    );
    if du == 0.0 && dv == 0.0 {
        return;
    }
    match &mut pcurve.geometry {
        PcurveGeometry::Nurbs(n) => {
            for p in &mut n.control_points {
                p.x += du;
                p.y += dv;
            }
        }
        PcurveGeometry::Elliptical { center, .. } => {
            center.0 += du;
            center.1 += dv;
        }
        PcurveGeometry::Linear { start, end } => {
            start.0 += du;
            start.1 += dv;
            end.0 += du;
            end.1 += dv;
        }
        _ => {}
    }
}

/// `P'(t) = P(t_start + t_end - t)` over the same parameter range.
fn reversed_pcurve(p: &Pcurve) -> Pcurve {
    let (a, b) = (p.t_start, p.t_end);
    let geometry = match &p.geometry {
        PcurveGeometry::Nurbs(n) => {
            let (k0, k1) = (n.knots[0], n.knots[n.knots.len() - 1]);
            let mut knots: Vec<f64> = n.knots.iter().map(|k| k0 + k1 - k).collect();
            knots.reverse();
            let mut control_points = n.control_points.clone();
            control_points.reverse();
            let mut weights = n.weights.clone();
            weights.reverse();
            PcurveGeometry::Nurbs(NurbsCurve {
                degree: n.degree,
                knots,
                control_points,
                weights,
            })
        }
        PcurveGeometry::Elliptical {
            center,
            u_radius,
            v_radius,
            rotation,
            start_angle,
            end_angle,
        } => PcurveGeometry::Elliptical {
            center: *center,
            u_radius: *u_radius,
            v_radius: *v_radius,
            rotation: *rotation,
            start_angle: *end_angle,
            end_angle: *start_angle,
        },
        PcurveGeometry::Linear { start, end } => PcurveGeometry::Linear {
            start: *end,
            end: *start,
        },
        other => other.clone(),
    };
    Pcurve {
        geometry,
        t_start: a,
        t_end: b,
    }
}

/// `surface(pcurve(t)) == curve(t)` (mirrored for a reversed pcurve)
/// at every verification sample.
fn verified(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    curve: &CurveGeometry,
    range: (f64, f64),
    reversed: bool,
    gate: f64,
) -> bool {
    (0..VERIFY_SAMPLES).all(|i| {
        let t = range.0 + (range.1 - range.0) * i as f64 / (VERIFY_SAMPLES - 1) as f64;
        let uv = pcurve.point_at(t);
        let on_surface = surface.point_at(uv.x, uv.y);
        let mirrored = if reversed { range.0 + range.1 - t } else { t };
        on_surface.distance_to(&curve.point_at(mirrored)) <= gate
    })
}
