//! Late-pipeline pass that clamps every half-edge's pcurve control
//! points to its surface's UV bounds for surfaces with bounded
//! parameter ranges (Gregory, NURBS, sphere, torus).
//!
//! Why a late pass: many cap-internal half-edges acquire pcurves
//! AFTER the cap surface is assigned — `attach_fillet_pcurves` runs
//! during fillet ribbon construction and inherits pcurves authored
//! against the parent (planar) face's UV space. By the time corner
//! blending runs `reproject_loop_pcurves_onto_gregory`, those
//! half-edges have no 3-D curve and no pcurve to fall back on, so
//! they're skipped — but later they get one, in the wrong UV box.
//!
//! This pass walks every half-edge with a pcurve on a bounded-UV
//! surface, samples the pcurve, and if any sample falls outside the
//! surface's parameter range it clamps every control point of the
//! pcurve into the box. Curve geometry may distort slightly near the
//! boundary but evaluations always land inside the valid domain —
//! the alternative is unreliable Gregory extrapolation that produces
//! 3-D points many model-units away from the actual edge.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::BRep;

/// Walk every half-edge with a pcurve on a bounded surface. If any
/// sample falls outside the surface's UV box, clamp every control
/// point to the box. Returns the number of pcurves clamped.
pub(crate) fn clamp_bounded_pcurves(brep: &mut BRep) -> usize {
    let n_samples = 9;
    let handles: Vec<popkern_core::topology::HalfEdgeHandle> = brep.half_edges.keys().collect();
    let mut clamped = 0;
    let mut total = 0;
    let mut with_pcurve = 0;
    let mut bounded_seen = 0;
    for heh in handles {
        total += 1;
        let Some(he_data) = brep.half_edges.get(heh) else {
            continue;
        };
        let Some(pcurve) = he_data.pcurve.as_ref() else {
            continue;
        };
        with_pcurve += 1;
        let Some(loop_data) = brep.loops.get(he_data.loop_) else {
            continue;
        };
        let face = loop_data.face;
        let Some(surf) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        // Only bounded surfaces — planes/cylinders use any-real (u, v)
        // and clamping doesn't help.
        let bounded = matches!(
            &surf.geometry,
            SurfaceGeometry::Nurbs(_)
                | SurfaceGeometry::Sphere { .. }
                | SurfaceGeometry::Torus { .. }
                | SurfaceGeometry::Gregory(_)
        );
        if !bounded {
            continue;
        }
        bounded_seen += 1;
        let ((u_lo, u_hi), (v_lo, v_hi)) = surf.geometry.parameter_range();
        if !u_lo.is_finite() || !u_hi.is_finite() || !v_lo.is_finite() || !v_hi.is_finite() {
            continue;
        }
        // Sphere u is 2π-periodic (atan2 wrap): a pcurve crossing
        // u = 0/2π unwraps to CPs outside [0, 2π], but `point_at`
        // uses cos(u)/sin(u) which evaluate correctly modulo 2π.
        // Clamping such CPs destroys the unwrap and produces 3-D
        // points up to a sphere-diameter away. Only clamp v on
        // Sphere; treat Cylinder similarly (its u is also 2π-
        // periodic, but Cylinder isn't bounded so it doesn't reach
        // here — handled by the `bounded` filter above).
        let u_periodic = matches!(&surf.geometry, SurfaceGeometry::Sphere { .. });
        // Check CONTROL POINTS directly, not curve evaluations: a
        // NURBS curve stays inside the convex hull of its CPs, so
        // any CP outside the box implies the curve can leave the
        // box (and *will* leave it for clamped knots — the first/
        // last CPs are exact endpoints). Sampling could miss
        // narrow excursions between samples, but CPs are the
        // ground truth.
        // Analytic Linear pcurves are in-box by construction (G4) —
        // EXCEPT for sub-tolerance endpoint drift from vertex/curve
        // repair passes (measured: pyramid4_zoo mixed served a
        // Linear endpoint at v = -1.05e-4 on a Gregory cap, and the
        // escape tessellated as a 3.1e-4 non-manifold seam). Clamp
        // Linear endpoints into the box when the drift is small
        // (<= 1e-3 of the range) — larger escapes are authoring
        // bugs that clamping would only mask.
        if let popkern_core::geometry::pcurve::PcurveGeometry::Linear { start, end } =
            &pcurve.geometry
        {
            let band_u = ((u_hi - u_lo).abs() * 1e-3).max(1e-9);
            let band_v = ((v_hi - v_lo).abs() * 1e-3).max(1e-9);
            let mut new_start = *start;
            let mut new_end = *end;
            let mut fix = false;
            for p in [&mut new_start, &mut new_end] {
                if !u_periodic {
                    if p.0 < u_lo && p.0 >= u_lo - band_u {
                        p.0 = u_lo;
                        fix = true;
                    } else if p.0 > u_hi && p.0 <= u_hi + band_u {
                        p.0 = u_hi;
                        fix = true;
                    }
                }
                if p.1 < v_lo && p.1 >= v_lo - band_v {
                    p.1 = v_lo;
                    fix = true;
                } else if p.1 > v_hi && p.1 <= v_hi + band_v {
                    p.1 = v_hi;
                    fix = true;
                }
            }
            if fix {
                let before = brep
                    .measure_half_edge_pcurve_deviation(heh)
                    .map(|r| r.max_deviation);
                let original = pcurve.clone();
                if let Some(he_mut) = brep.half_edges.get_mut(heh)
                    && let Some(pc) = he_mut.pcurve.as_mut()
                {
                    pc.geometry = popkern_core::geometry::pcurve::PcurveGeometry::Linear {
                        start: new_start,
                        end: new_end,
                    };
                }
                let after = brep
                    .measure_half_edge_pcurve_deviation(heh)
                    .map(|r| r.max_deviation);
                if std::env::var("POPKERN_CLAMP_TRACE").is_ok() {
                    eprintln!(
                        "  [linear_clamp] he {heh:?} face {face:?} dev {before:?} -> {after:?}"
                    );
                }
                match (before, after) {
                    (Some(b), Some(a)) if a <= b + 1e-9 => clamped += 1,
                    _ => {
                        if let Some(he_mut) = brep.half_edges.get_mut(heh) {
                            he_mut.pcurve = Some(original);
                        }
                    }
                }
            }
            continue;
        }
        let Some(nurbs_view) = pcurve.as_nurbs() else {
            continue;
        };
        let cps_outside = nurbs_view.control_points.iter().any(|p| {
            let u_out = if u_periodic {
                false
            } else {
                p.x < u_lo - 1e-9 || p.x > u_hi + 1e-9
            };
            u_out || p.y < v_lo - 1e-9 || p.y > v_hi + 1e-9
        });
        if std::env::var("POPKERN_CLAMP_TRACE").is_ok() {
            let kind = match &surf.geometry {
                SurfaceGeometry::Gregory(_) => "Gregory",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                _ => "other",
            };
            let cp_summary: Vec<String> = nurbs_view
                .control_points
                .iter()
                .map(|p| format!("({:.2},{:.2})", p.x, p.y))
                .collect();
            eprintln!(
                "  [gregory_clamp_inspect] he {heh:?} face {face:?} ({kind}) cps_outside={cps_outside} cps={}",
                cp_summary.join(",")
            );
        }
        if !cps_outside {
            continue;
        }
        let before_deviation = brep
            .measure_half_edge_pcurve_deviation(heh)
            .map(|report| report.max_deviation);
        let original_pcurve = pcurve.clone();
        let _ = n_samples;
        // Clamp every CP into the box. `Some(...)`-mut access via
        // brep.half_edges to update.
        let Some(he_mut) = brep.half_edges.get_mut(heh) else {
            continue;
        };
        if let Some(pc) = he_mut.pcurve.as_mut()
            && let Some(n) = pc.as_nurbs_mut()
        {
            for cp in &mut n.control_points {
                if !u_periodic {
                    cp.x = cp.x.clamp(u_lo, u_hi);
                }
                cp.y = cp.y.clamp(v_lo, v_hi);
            }
        }
        if clamp_worsened_edge_agreement(brep, heh, before_deviation) {
            if let Some(he_mut) = brep.half_edges.get_mut(heh) {
                he_mut.pcurve = Some(original_pcurve);
            }
            continue;
        }
        clamped += 1;
    }
    if std::env::var("POPKERN_CLAMP_TRACE").is_ok() {
        eprintln!(
            "  [gregory_pcurve_clamp] total_he={total} with_pcurve={with_pcurve} bounded_seen={bounded_seen} clamped={clamped}"
        );
    }
    clamped
}

fn clamp_worsened_edge_agreement(
    brep: &BRep,
    heh: popkern_core::topology::HalfEdgeHandle,
    before_deviation: Option<f64>,
) -> bool {
    let Some(before) = before_deviation else {
        return false;
    };
    let Some(after) = brep
        .measure_half_edge_pcurve_deviation(heh)
        .map(|report| report.max_deviation)
    else {
        return true;
    };
    if !after.is_finite() {
        return true;
    }
    after > before + brep.tolerance.linear().max(1.0e-9)
}
