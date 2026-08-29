//! EDGE_CURVE `same_sense` derivation for STEP export.
//!
//! Split from `export_periodic_branches.rs` to honour the structural
//! line ratchet. Uses that module's chain-direction normalization
//! for closed rims.

use popkern_core::topology::BRep;

use super::export_periodic_branches::{
    PeriodicNormalization, export_pcurve_for_use, pcurve_opposes_edge_curve,
};

/// EDGE_CURVE same_sense: whether the edge's topological direction
/// — defined by its CANONICAL half-edge, the one every loop writes
/// as the FORWARD ORIENTED_EDGE — agrees with the basis curve's
/// parameterization. STEP consumers compute each use's traversal as
/// curve_direction × same_sense × oriented_edge_sense; deriving the
/// flag from anything other than the canonical half-edge (the old
/// code used the trim ordering t_end ≥ t_start) breaks that product
/// for one side of every edge whose canonical use runs against the
/// curve. OpenCASCADE reprojection-heals the damage on simple faces
/// but not on periodic ones: the flagship seam-rim circles are CCW
/// curves whose canonical use is the CW plate-hole side, and the
/// wrongly-true flag reversed the torus-side traversal by a full
/// period (self-intersecting wire, the FreeCAD sliver).
pub(super) fn edge_same_sense(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
    normalization: &PeriodicNormalization,
) -> bool {
    use popkern_core::geometry::CurveEval;
    let Some(edge_data) = brep.edges.get(edge) else {
        return true;
    };
    let Some(curve) = edge_data.curve.as_ref() else {
        return true;
    };
    let canonical = edge_data.half_edge;
    let Some(he) = brep.half_edges.get(canonical) else {
        return true;
    };
    let ceval: &dyn CurveEval = &curve.geometry;
    let start = ceval.point_at(curve.t_start);
    let end = ceval.point_at(curve.t_end);
    // same_sense relates the edge direction to the curve's UNDERLYING
    // parameterization; a decreasing trim (t_end < t_start) runs the
    // edge against it, so every direction test XORs the trim order.
    let trim_increasing = curve.t_end >= curve.t_start;
    if start.distance_to(&end) > 1e-9 {
        // Open edge: the canonical half-edge's origin sits at the
        // trim end it departs from.
        let Some(origin) = brep.vertices.get(he.origin).map(|v| v.position) else {
            return true;
        };
        let origin_at_trim_start = origin.distance_to(&start) <= origin.distance_to(&end);
        return origin_at_trim_start == trim_increasing;
    }
    // Closed rim: endpoints coincide, so the canonical direction
    // must come from a loop-level signal. Prefer the chain direction
    // the periodic normalization resolved by loop closure (the
    // manifold-consistent traversal); fall back to the stored
    // pcurve's orientation (traversal-aligned by contract, not
    // enforced — PK-100) for uses outside periodic loops.
    let face_handle = match brep.loops.get(he.loop_) {
        Some(l) => l.face,
        None => return true,
    };
    let Some(surface) = brep.faces.get(face_handle).and_then(|f| f.surface.as_ref()) else {
        return true;
    };
    let Some(pcurve) = export_pcurve_for_use(brep, face_handle, canonical, &surface.geometry)
    else {
        return true;
    };
    let sense_trace = popkern_core::trace::trace_enabled("step-sense");
    if let Some(&chain_forward) = normalization.directions.get(&(face_handle, canonical)) {
        let opposes = pcurve_opposes_edge_curve(brep, canonical, &surface.geometry, &pcurve);
        let along_trim = chain_forward != opposes;
        if sense_trace {
            eprintln!(
                "[sense] {edge:?} closed-rim CHAIN branch: face {face_handle:?} chain_forward {chain_forward} opposes {opposes} trim_inc {trim_increasing} -> {}",
                along_trim == trim_increasing
            );
        }
        return along_trim == trim_increasing;
    }
    // No chain info for the canonical use (its face is not periodic
    // — e.g. a plane's single-edge hole loop). The TWIN use may be
    // chain-resolved; the canonical direction is its opposite.
    let twin = he.twin;
    if let Some(twin_he) = brep.half_edges.get(twin)
        && let Some(twin_loop) = brep.loops.get(twin_he.loop_)
        && let Some(&twin_forward) = normalization.directions.get(&(twin_loop.face, twin))
        && let Some(twin_surface) = brep
            .faces
            .get(twin_loop.face)
            .and_then(|f| f.surface.as_ref())
        && let Some(twin_pcurve) =
            export_pcurve_for_use(brep, twin_loop.face, twin, &twin_surface.geometry)
    {
        let twin_opposes =
            pcurve_opposes_edge_curve(brep, twin, &twin_surface.geometry, &twin_pcurve);
        let twin_along_trim = twin_forward != twin_opposes;
        if sense_trace {
            eprintln!(
                "[sense] {edge:?} closed-rim TWIN branch: twin_forward {twin_forward} twin_opposes {twin_opposes} trim_inc {trim_increasing} -> {}",
                !twin_along_trim == trim_increasing
            );
        }
        return !twin_along_trim == trim_increasing;
    }
    // Last resort: the stored pcurve's orientation (traversal-aligned
    // by contract, not enforced — PK-100).
    let opposes = pcurve_opposes_edge_curve(brep, canonical, &surface.geometry, &pcurve);
    if sense_trace {
        eprintln!(
            "[sense] {edge:?} closed-rim LAST-RESORT branch: face {face_handle:?} opposes {opposes} trim_inc {trim_increasing} -> {}",
            !opposes == trim_increasing
        );
    }
    !opposes == trim_increasing
}
