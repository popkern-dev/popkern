//! Import precision tracking (split from `import.rs`, 500-line cap):
//! measured edge deviations become edge tolerances, curve-end gaps
//! within the file's declared precision become vertex tolerances
//! (PK-639), and the body's precision envelope widens to cover them.

use crate::step::diagnostics::{DiagnosticLog, Phase};

/// Measure edge deviations on an imported BRep and set entity tolerances.
///
/// For each edge with a curve, measures the deviation between the 3D curve
/// and the surface evaluations via pcurves. If the deviation exceeds
/// a small threshold, records it as the edge's tracked imprecision.
pub(crate) fn measure_import_deviations(
    brep: &mut popkern_core::topology::BRep,
    record_above: f64,
    end_gap_allowance: f64,
    log: &mut DiagnosticLog,
) {
    let edge_handles: Vec<_> = brep.edges.keys().collect();
    // A measured deviation is a fact about the file, recorded whatever
    // the tolerance policy: recording only those above the body's
    // tolerance dropped every deviation a looser file uncertainty
    // admitted, and seam certification then read them as open gaps
    // (measured: OCC's rewrite of offset_bracket declares 2e-3; seven
    // seams certified G0 at the default and "open" at 2e-3).
    let linear_tol = record_above;
    let mut edge_count = 0u32;
    let mut max_tracked_deviation = 0.0_f64;

    // Curve ends off their vertices within the file's declared
    // precision record that gap on the VERTICES (a curve end within a
    // vertex's tolerance is that vertex); a larger gap is a defect (a
    // corrupted coordinate) the validator must see.
    for &eh in &edge_handles {
        record_end_gap(brep, eh, end_gap_allowance);
    }
    for eh in edge_handles {
        // An edge's tolerance covers its own deviation from its faces
        // and, by the kernel's containment rule (vertex <= edge), the
        // tolerances of its end vertices.
        let ends = brep
            .half_edges
            .get(brep.edges[eh].half_edge)
            .and_then(|he| Some((he.origin, brep.half_edges.get(he.twin)?.origin)));
        let vertex_tol = ends.map_or(0.0, |(a, b)| {
            brep.vertex_tolerance(a).max(brep.vertex_tolerance(b))
        });
        let dev = brep.measure_edge_deviation(eh).max(vertex_tol);
        if dev > linear_tol && dev > brep.edge_tolerance(eh) {
            brep.set_entity_tolerance(popkern_core::precision::EntityHandle::Edge(eh), dev);
            max_tracked_deviation = max_tracked_deviation.max(dev);
            edge_count += 1;
        }
    }

    widen_import_precision_envelope(brep, max_tracked_deviation);

    if edge_count > 0 {
        log.warn(
            Phase::Validate,
            None,
            format!(
                "{edge_count} edge(s) with measurable deviation tracked; import precision envelope {:.2e}",
                brep.max_imprecision
            ),
        );
    }
}

/// Record an edge curve's end gaps (curve ends against its two
/// vertices, as an unordered pair) on the vertices when within
/// `allowance`.
fn record_end_gap(
    brep: &mut popkern_core::topology::BRep,
    eh: popkern_core::topology::EdgeHandle,
    allowance: f64,
) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::precision::EntityHandle;
    let ends = || -> Option<_> {
        let edge = brep.edges.get(eh)?;
        let curve = edge.curve.as_ref()?;
        let he = brep.half_edges.get(edge.half_edge)?;
        let twin = brep.half_edges.get(he.twin)?;
        let (va, vb) = (he.origin, twin.origin);
        let (a, b) = (
            brep.vertices.get(va)?.position,
            brep.vertices.get(vb)?.position,
        );
        let (s, e) = (
            curve.geometry.point_at(curve.t_start),
            curve.geometry.point_at(curve.t_end),
        );
        let natural = (a.distance_to(&s), b.distance_to(&e));
        let swapped = (a.distance_to(&e), b.distance_to(&s));
        let gaps = if natural.0.max(natural.1) <= swapped.0.max(swapped.1) {
            natural
        } else {
            swapped
        };
        Some((va, vb, gaps))
    };
    let Some((va, vb, (ga, gb))) = ends() else {
        return;
    };
    for (v, gap) in [(va, ga), (vb, gb)] {
        if gap > brep.vertex_tolerance(v) && gap <= allowance {
            brep.set_entity_tolerance(EntityHandle::Vertex(v), gap);
        }
    }
}

fn widen_import_precision_envelope(
    brep: &mut popkern_core::topology::BRep,
    max_tracked_deviation: f64,
) {
    if max_tracked_deviation <= brep.max_imprecision {
        return;
    }
    brep.max_imprecision = (max_tracked_deviation * 1.01)
        .max(max_tracked_deviation + brep.tolerance.numerical_floor());
}
