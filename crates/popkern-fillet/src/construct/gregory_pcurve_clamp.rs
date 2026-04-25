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
        let face = brep.loops[he_data.loop_].face;
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
        // Check CONTROL POINTS directly, not curve evaluations: a
        // NURBS curve stays inside the convex hull of its CPs, so
        // any CP outside the box implies the curve can leave the
        // box (and *will* leave it for clamped knots — the first/
        // last CPs are exact endpoints). Sampling could miss
        // narrow excursions between samples, but CPs are the
        // ground truth.
        let cps_outside = pcurve.control_points.iter().any(|p| {
            p.x < u_lo - 1e-9 || p.x > u_hi + 1e-9 || p.y < v_lo - 1e-9 || p.y > v_hi + 1e-9
        });
        if std::env::var("POPKERN_CLAMP_TRACE").is_ok() {
            let kind = match &surf.geometry {
                SurfaceGeometry::Gregory(_) => "Gregory",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                _ => "other",
            };
            let cp_summary: Vec<String> = pcurve
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
        let _ = n_samples;
        // Clamp every CP into the box. `Some(...)`-mut access via
        // brep.half_edges to update.
        let Some(he_mut) = brep.half_edges.get_mut(heh) else {
            continue;
        };
        if let Some(pc) = he_mut.pcurve.as_mut() {
            for cp in &mut pc.control_points {
                cp.x = cp.x.clamp(u_lo, u_hi);
                cp.y = cp.y.clamp(v_lo, v_hi);
            }
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
