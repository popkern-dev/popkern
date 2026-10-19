//! Trim pcurves of a closed fillet's rails on their analytic supports
//! (PK-292): the exact plane arm or the adaptive chart fit, angle-
//! matched to the rail. A NURBS face's pcurves are the projection
//! recompute's (PK-891).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle};

/// Author the pcurve of a closed-fillet trim half-edge for its exact
/// `Arc` rail (step 13b of the closed fillet) on an analytic face: the
/// exact plane arm, else the adaptive chart fit (`fit_chart_pcurve`),
/// both parameterised by the SAME angles as the arc and verified.
/// Without angle-matched pcurves the STEP exporter's pcurve synthesis
/// disagrees with the exact Arc and the wall face imports broken in
/// OpenCASCADE (PK-292). Leaves the half-edge untouched when every
/// author declines.
pub(crate) fn author_trim_pcurve(brep: &mut BRep, he: HalfEdgeHandle, tolerance: &Tolerance) {
    let (edge_h, face_h) = {
        let d = &brep.half_edges[he];
        (d.edge, brep.loops[d.loop_].face)
    };
    let Some(curve) = brep.edges[edge_h].curve.clone() else {
        return;
    };
    let Some(sof) = brep.faces[face_h].surface.clone() else {
        return;
    };
    // A NURBS face's trim pcurves are left to the projection recompute
    // (`BRep::recompute_face_pcurves`), whose branch chaining runs each
    // closed rail one period round a periodic strip (PK-891) and whose
    // pcurves the seam certifier certifies (the G2 rim fixture). A
    // boundary-iso search authored them here until PK-891: on the twin
    // side of a closed edge it ran the curve rather than the half-edge,
    // so the stale repair rebuilt every one anyway.
    if matches!(sof.geometry, SurfaceGeometry::Nurbs(_)) {
        return;
    }
    let range = (curve.t_start, curve.t_end);
    if let Some(pc) = popkern_core::geometry::chart_pcurve::analytic_chart_pcurve(
        &sof.geometry,
        &curve.geometry,
        range,
    ) {
        brep.half_edges[he].pcurve = Some(pc);
        return;
    }
    if let Some(pc) = popkern_core::geometry::chart_pcurve::fit_chart_pcurve(
        &sof.geometry,
        &curve.geometry,
        range,
        tolerance,
    ) {
        brep.half_edges[he].pcurve = Some(pc.into());
    }
}
