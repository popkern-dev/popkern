//! Corner-resolution gate (PK-516) and the boundary-ring fallback's
//! cap sense (PK-517), pinned on the pentagonal-pyramid apex.
//!
//! The setback vertex blend serves this apex in production, so the
//! fallback ladder behind it is reached here through the thread-scoped
//! `with_setback_forced_to_decline` test seam — the stand-in for "any
//! future setback decline".

use popkern_blend::error::{BlendError, BlendFailure};
use popkern_blend::{
    ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, VertexMode, fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

fn pentagon_apex() -> (BRep, VertexHandle, Vec<EdgeHandle>) {
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, Tolerance::default()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges = pyramid
        .edges
        .iter()
        .filter_map(|(eh, edge)| {
            let he = &pyramid.half_edges[edge.half_edge];
            let twin_origin = pyramid.half_edges[he.twin].origin;
            (he.origin == apex || twin_origin == apex).then_some(eh)
        })
        .collect();
    (pyramid, apex, apex_edges)
}

fn apex_fillet(vertex_mode: VertexMode) -> (BRep, VertexHandle, Result<(), BlendError>) {
    let (pyramid, apex, edges) = pentagon_apex();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.16),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        vertex_mode,
        ..Default::default()
    };
    let result = popkern_blend::dispatch_overrides::with_setback_forced_to_decline(|| {
        fillet_edges(&pyramid, &specs, &options, Tolerance::default(), None).map(|_| ())
    });
    (pyramid, apex, result)
}

/// The decline must be the typed corner failure at the apex, and the
/// operation transactional: the partial result is the untouched input.
fn corner_decline(input: &BRep, apex: VertexHandle, result: Result<(), BlendError>) -> String {
    let Err(BlendError::PartialFailure {
        partial_result,
        failures,
        ..
    }) = result
    else {
        panic!("an unresolvable corner must decline the fillet, got {result:?}");
    };
    assert_eq!(
        partial_result.faces.len(),
        input.faces.len(),
        "a declined fillet hands back the untouched input"
    );
    let [
        (
            _,
            BlendFailure::VertexBlendFailed {
                vertex,
                reason,
                location,
            },
        ),
    ] = failures.as_slice()
    else {
        panic!("expected exactly one VertexBlendFailed, got {failures:?}");
    };
    assert_eq!(*vertex, apex);
    let at = location.expect("the decline carries the corner position");
    assert!(at.distance_to(&input.vertices[apex].position) < 1e-9);
    assert!(
        reason.contains("curve-surface consistency"),
        "the decline names the measured defect: {reason}"
    );
    reason.clone()
}

/// PK-516: before the gate, this ladder ended in a mitre whose ten rim
/// edges sat 1.87e-2 off their strips, returned as `Ok` because the
/// fast `validate()` excludes check 9. Now it declines, and the reason
/// carries every cap route's verdict.
///
/// PK-517 rides the same ladder: the boundary-ring rung builds the
/// cap on a pushed ring. Before the ring took its face-loop sense
/// from the solid, and before PnS aligned n >= 5 sectors to the loop
/// and read the tensor border inward, the PnS cap landed every
/// sector on the wrong sub-face (rim g0 2.6e-1, pole 0.9528 BELOW
/// the 0.9603 rim) and the N-quad Gregory cap was mirrored (rim G1
/// normal deficit 0.999). Both still decline on this oversized
/// pushed ring — PnS's quadratic border misses the rim arcs by
/// 2.3e-2 and the N-quad's centroid pole cannot turn tangent to the
/// strips (0.58) — but by honest construction limits, not by an
/// inverted cap.
#[test]
fn smooth_fallback_declines_instead_of_serving_an_off_surface_mitre() {
    let (input, apex, result) = apex_fillet(VertexMode::Smooth);
    let reason = corner_decline(&input, apex, result);
    assert!(
        reason.contains("all-or-nothing rollback"),
        "the N=5 corner rolled back: {reason}"
    );
    let measured = |stage: &str, metric: &str| -> f64 {
        let verdict = reason
            .split("; ")
            .find(|v| v.contains(stage))
            .unwrap_or_else(|| panic!("no {stage} verdict in: {reason}"));
        let after = &verdict[verdict.find(metric).expect("metric") + metric.len()..];
        after
            .split_whitespace()
            .next()
            .and_then(|v| v.parse().ok())
            .expect("numeric metric")
    };
    let pns_rim_g0 = measured("pns.dispatch.cap_support_acceptance", "cap_outer_g0");
    assert!(
        pns_rim_g0 < 5.0e-2,
        "PnS sectors must land on their own sub-faces (was 2.6e-1): {pns_rim_g0:.3e}"
    );
    let gregory_rim_g1 = measured(
        "rational_gregory.dispatch.cap_support_acceptance",
        "cap_outer_g1",
    );
    assert!(
        gregory_rim_g1 < 0.9,
        "the N-quad cap must not be mirrored (was 0.999): {gregory_rim_g1:.3e}"
    );
}

/// PK-516 + PK-522: an explicit mitre request at this G2 apex keeps
/// its off-surface residue. The mitre-corner settlement trims only
/// circular-cylinder strips, and curvature-continuous G2 strips are
/// not cylinders. The request declines with both causes named.
#[test]
fn requested_mitre_declines_when_its_corner_residue_is_off_surface() {
    let (input, apex, result) = apex_fillet(VertexMode::Mitred);
    let reason = corner_decline(&input, apex, result);
    assert!(
        reason.contains("mitred corner requested; mitre declined: a strip is not a cylinder"),
        "the decline names the cause: {reason}"
    );
}
