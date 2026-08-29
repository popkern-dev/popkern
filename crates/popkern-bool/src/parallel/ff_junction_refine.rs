//! Junction-endpoint boundary refinement (PK-313).
//!
//! Per-pair clips evaluate FITTED record curves at face-boundary
//! crossings, so each pair's junction endpoint carries the fit's
//! lateral slop — measured up to ~1e-3 at grazing junction clusters,
//! where coincident boundary events split apart, the connector arc
//! collapses below tolerance, and assembly declines on a half-claimed
//! strand (the outer-grazing torus micro-lens). Refine every record
//! endpoint that sits near a face iso-boundary onto the TRUE
//! surface×surface intersection ∩ boundary point via the pinned
//! Newton (`refine_boundary_crossing`): agreeing pairs then land on
//! the same junction bit-exactly, the way seam-refined junctions
//! already do.
//!
//! Honesty bounds: a refinement is accepted only when the Newton
//! converges (a curve tangent to the boundary goes singular and keeps
//! its unrefined endpoint) and the endpoint moves less than a quarter
//! of the record's chord. Endpoints far from every iso-boundary (free
//! curve ends inside a face) are left untouched.

use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::geometry::intersection::marching::{BoundaryPin, refine_boundary_crossing};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::intersect::IntersectionRecord;

/// Chart-distance band within which an endpoint counts as sitting on
/// a face iso-boundary and earns refinement. Generous relative to the
/// measured slop class (~1e-3–3e-3): refining an already-exact
/// endpoint is a converged no-op, so the band's only cost is a Newton
/// solve per qualifying endpoint.
const PIN_BAND: f64 = 0.05;

/// Nearest-image distance of `x` to `value` under an optional period.
fn axis_gap(x: f64, value: f64, period: Option<f64>) -> f64 {
    let d = x - value;
    match period {
        Some(p) if p.is_finite() && p > 0.0 => d - (d / p).round() * p,
        _ => d,
    }
}

/// Refine the endpoints of freshly clipped per-pair records onto their
/// nearest face iso-boundaries. See the module docs for the defect
/// class and the honesty bounds.
#[allow(clippy::too_many_arguments)]
pub(super) fn refine_record_junction_endpoints(
    records: &mut [IntersectionRecord],
    a: &BRep,
    face_a: FaceHandle,
    surf_a: &SurfaceGeometry,
    b: &BRep,
    face_b: FaceHandle,
    surf_b: &SurfaceGeometry,
    tolerance: &Tolerance,
) {
    if records.is_empty() {
        return;
    }
    let cell_a = crate::clip::chart_cell::face_chart_cell(a, face_a, tolerance);
    let cell_b = crate::clip::chart_cell::face_chart_cell(b, face_b, tolerance);
    if cell_a.is_none() && cell_b.is_none() {
        return;
    }
    let (pu_a, pv_a) = surf_a.periodicity();
    let (pu_b, pv_b) = surf_b.periodicity();
    let trace = popkern_core::trace::trace_enabled("boolean");

    for record in records.iter_mut() {
        let (t0, t1) = record.curve.curve_3d.parameter_range();
        let chord = record
            .curve
            .curve_3d
            .point_at(t0)
            .distance_to(&record.curve.curve_3d.point_at(t1));
        let move_bound = (0.25 * chord).max(tolerance.linear() * 100.0);
        for (end, t) in [(false, t0), (true, t1)] {
            let p = record.curve.curve_3d.point_at(t);
            let uv1 = record.curve.curve_uv1.point_at(t);
            let uv2 = record.curve.curve_uv2.point_at(t);
            // Candidate pins: nearest iso-edge of each chart cell.
            let mut best: Option<(f64, BoundaryPin)> = None;
            let mut consider = |on_surface_one: bool,
                                axis_u: bool,
                                coord: f64,
                                value: f64,
                                period: Option<f64>| {
                let gap = axis_gap(coord, value, period).abs();
                if gap < PIN_BAND && best.as_ref().is_none_or(|(bg, _)| gap < *bg) {
                    // Pin value on the seed's own period branch.
                    let value_here = coord - axis_gap(coord, value, period);
                    best = Some((
                        gap,
                        BoundaryPin {
                            on_surface_one,
                            axis_u,
                            value: value_here,
                        },
                    ));
                }
            };
            if let Some(cell) = &cell_a {
                consider(true, true, uv1.x, cell.u_range.0, pu_a);
                consider(true, true, uv1.x, cell.u_range.1, pu_a);
                consider(true, false, uv1.y, cell.v_range.0, pv_a);
                consider(true, false, uv1.y, cell.v_range.1, pv_a);
            }
            if let Some(cell) = &cell_b {
                consider(false, true, uv2.x, cell.u_range.0, pu_b);
                consider(false, true, uv2.x, cell.u_range.1, pu_b);
                consider(false, false, uv2.y, cell.v_range.0, pv_b);
                consider(false, false, uv2.y, cell.v_range.1, pv_b);
            }
            let Some((gap, pin)) = best else { continue };
            if gap <= 1.0e-10 {
                continue; // already exact
            }
            let Some((u1r, v1r, u2r, v2r, refined)) = refine_boundary_crossing(
                surf_a, uv1.x, uv1.y, surf_b, uv2.x, uv2.y, pin, tolerance,
            ) else {
                continue; // tangent/singular — keep unrefined endpoint
            };
            let moved = refined.distance_to(&p);
            if moved > move_bound {
                continue;
            }
            // Write back on each chart's original period branch.
            let branch = |x: f64, reference: f64, period: Option<f64>| match period {
                Some(p) if p.is_finite() && p > 0.0 => x - ((x - reference) / p).round() * p,
                _ => x,
            };
            let u1w = branch(u1r, uv1.x, pu_a);
            let v1w = branch(v1r, uv1.y, pv_a);
            let u2w = branch(u2r, uv2.x, pu_b);
            let v2w = branch(v2r, uv2.y, pv_b);
            if trace {
                eprintln!(
                    "  [junction-refine] ({face_a:?},{face_b:?}) {} gap {gap:.2e} moved {moved:.2e} -> ({:.6},{:.6},{:.6})",
                    if end { "end" } else { "start" },
                    refined.x,
                    refined.y,
                    refined.z
                );
            }
            let set =
                |c: &mut popkern_core::geometry::NurbsCurve, x: f64, y: f64, z: Option<f64>| {
                    let cp = if end {
                        c.control_points.last_mut()
                    } else {
                        c.control_points.first_mut()
                    };
                    if let Some(cp) = cp {
                        cp.x = x;
                        cp.y = y;
                        if let Some(z) = z {
                            cp.z = z;
                        }
                    }
                };
            set(
                &mut record.curve.curve_3d,
                refined.x,
                refined.y,
                Some(refined.z),
            );
            set(&mut record.curve.curve_uv1, u1w, v1w, None);
            set(&mut record.curve.curve_uv2, u2w, v2w, None);
        }
    }
}
