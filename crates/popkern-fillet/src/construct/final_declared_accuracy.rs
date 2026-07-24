//! Final declared-accuracy pass (PK-177).
//!
//! After every geometry author has run, an edge can carry a small
//! honest residual to its adjacent surfaces that no refit can
//! reduce (surface sampling density, independently fitted strip
//! halves, settled endpoints). Recording that measured residual as
//! the edge's tolerance is the honest accuracy declaration — the
//! serving gate and the STEP export both judge coherence at the
//! declared accuracy. Only SUB-PRECISION residuals record (< 5e-5,
//! margin capped at 9e-5 under the kernel's 1e-4 Precision
//! ceiling); anything larger is a real defect left exposed.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::BRep;

/// Measure each edge's worst curve-to-adjacent-surface distance
/// (continuation-seeded projection) and record sub-precision
/// residuals above the edge's current band as its tolerance.
/// Returns how many edges recorded.
pub(crate) fn record_final_declared_accuracy(brep: &mut BRep) -> usize {
    let model_band = brep.tolerance.linear().max(1e-9);
    let trace = std::env::var_os("PK_RECONCILE_TRACE").is_some();
    let mut recorded = 0usize;
    let mut measured = 0usize;
    let mut skipped_large = 0usize;
    let edges: Vec<_> = brep.edges.keys().collect();
    for edge in edges {
        let Some(ed) = brep.edges.get(edge) else {
            continue;
        };
        let Some(curve) = ed.curve.clone() else {
            continue;
        };
        let current = ed.tolerance.unwrap_or(0.0);
        let Some(he) = brep.half_edges.get(ed.half_edge) else {
            continue;
        };
        let Some(tw) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let faces = [
            brep.loops.get(he.loop_).map(|l| l.face),
            brep.loops.get(tw.loop_).map(|l| l.face),
        ];
        let mut worst = 0.0_f64;
        for face in faces.into_iter().flatten() {
            let Some(sof) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let mut prev_uv: Option<(f64, f64)> = None;
            for i in 0..=16 {
                let t = curve.t_start + (curve.t_end - curve.t_start) * (i as f64) / 16.0;
                let pt = curve.geometry.point_at(t);
                // Per-sample: best of the cold and continuation-
                // seeded projections (either alone can clamp to the
                // wrong chart region on extended strips).
                let cold =
                    popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &pt);
                let seeded = prev_uv.and_then(|(pu, pv)| {
                    popkern_core::geometry::uv_projection::project_point_to_uv_near(
                        &sof.geometry,
                        &pt,
                        pu,
                        pv,
                    )
                });
                let best = [cold, seeded]
                    .into_iter()
                    .flatten()
                    .map(|(u, v)| (pt.distance_to(&sof.geometry.point_at(u, v)), (u, v)))
                    .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                if let Some((d, uv)) = best {
                    prev_uv = Some(uv);
                    worst = worst.max(d);
                }
            }
        }
        measured += 1;
        if worst >= 5.0e-5 {
            skipped_large += 1;
        }
        // Record only genuinely sub-precision residuals above the
        // edge's current declared band.
        if worst > model_band.max(current)
            && worst < 5.0e-5
            && let Some(edge_mut) = brep.edges.get_mut(edge)
        {
            edge_mut.tolerance = Some(current.max((worst * 1.5).min(9.0e-5)));
            recorded += 1;
            if std::env::var_os("PK_RECONCILE_TRACE").is_some() {
                eprintln!(
                    "[declared-accuracy] {edge:?} residual {worst:.3e} recorded {:.3e}",
                    (worst * 1.5).min(9.0e-5)
                );
            }
        }
    }
    if trace {
        eprintln!(
            "[declared-accuracy] measured {measured} edge(s), recorded {recorded}, \
             above-ceiling {skipped_large}"
        );
    }
    recorded
}
