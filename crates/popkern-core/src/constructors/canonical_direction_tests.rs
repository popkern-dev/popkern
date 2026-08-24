//! Closed-edge canonical direction pins for constructors.

use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Closed-edge canonical convention (see `recompute_face_pcurves`): the
/// edge's `half_edge` must traverse the stored curve FORWARD. Verified
/// by comparing each canonical half's authored pcurve direction (mapped
/// through its face surface) against the curve tangent at t_start. A
/// violation is invisible on the constructor output (authored pcurves
/// are kept) but flips loop windings whenever pcurves are re-authored
/// (shell/offset assembly re-derives every pcurve).
#[test]
fn closed_edge_canonical_half_traverses_curve_forward() {
    use crate::geometry::curve::CurveEval;
    use crate::geometry::surface::SurfaceEval;

    for (name, brep) in [
        (
            "cylinder",
            BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap(),
        ),
        ("cone", BRep::cone(1.0, 2.0, Tolerance::default()).unwrap()),
    ] {
        for (eh, edge) in &brep.edges {
            let Some(curve) = &edge.curve else { continue };
            let start = curve.geometry.point_at(curve.t_start);
            let end = curve.geometry.point_at(curve.t_end);
            if start.distance_to(&end) > 1e-9 {
                continue; // open edge: direction disambiguated by endpoints
            }
            let he = &brep.half_edges[edge.half_edge];
            let Some(pc) = &he.pcurve else { continue };
            let face = brep.loops[he.loop_].face;
            let surface = &brep.faces[face].surface.as_ref().unwrap().geometry;
            // Map the canonical half's pcurve start tangent to 3D and
            // compare with the curve's forward tangent at t_start.
            let (p0, p1) = (
                pc.point_at(pc.t_start),
                pc.point_at(pc.t_start + (pc.t_end - pc.t_start) * 1e-3),
            );
            let s0 = surface.point_at(p0.x, p0.y);
            let s1 = surface.point_at(p1.x, p1.y);
            let he_dir = s1 - s0;
            let curve_dir = curve.geometry.tangent_at(curve.t_start);
            assert!(
                he_dir.dot(&curve_dir) > 0.0,
                "{name} edge {eh:?}: canonical half traverses its closed \
                 curve BACKWARD (convention: edge.half_edge = forward)"
            );
        }
    }
}
