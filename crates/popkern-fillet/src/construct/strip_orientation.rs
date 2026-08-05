//! Strip-face orientation honesty at commit (PK-238).
//!
//! The fillet pipeline authored strip surfaces with whatever natural
//! normal direction the construction produced and stamped
//! `orientation: true` unconditionally. On outward-oriented shells
//! (boolean/sweep output) that happened to be honest; on
//! inward-oriented shells (primitive-constructor bodies such as
//! `BRep::unit_cube`) the strip face became the ONLY
//! flag-inconsistent face in its shell. Every consumer that trusts
//! stored orientation flags — edge convexity classification,
//! effective-outward folding, rolling-ball side selection — then
//! broke on the NEXT operation that read the strip: the §12.3
//! blend-on-blend contact chain classified the strip-end arc
//! Concave (it is a 90° convex wedge) and solved the non-material
//! offset branch (PK-238: strip shipped on the y=1.02/R=0.08 torus
//! instead of the material 0.98/0.08).
//!
//! The reconcile rule is convention-independent and local: across a
//! tangency trim the strip and its support face are G1-tangent, and
//! one manifold shell requires their ORIENTED normals to agree
//! there. The strip's flag is therefore chosen so its flag-folded
//! natural normal matches the support face's flag-folded natural
//! normal at a trim-edge sample. Transverse neighbours (end caps,
//! corner bridges) are rejected by the |cos| witness floor; with no
//! tangent witness at all the authored flag is left unchanged.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Minimum |cos| between the strip normal and a support face's
/// oriented normal for the edge to count as a tangency witness.
/// Genuine trim contacts measure ~1 within construction tolerance;
/// transverse end-cap and bridge edges measure ~0.
const TANGENCY_WITNESS_FLOOR: f64 = 0.9;

/// Set `strip_face`'s orientation flag so its oriented normal agrees
/// with its support faces across the given tangency trim edges.
/// Returns the applied flag, or `None` when no edge produced a
/// usable tangency witness (flag left as authored).
pub(crate) fn reconcile_strip_orientation(
    brep: &mut BRep,
    strip_face: FaceHandle,
    trim_edges: &[EdgeHandle],
) -> Option<bool> {
    let strip_geometry = brep
        .faces
        .get(strip_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())?;

    let mut best: Option<(f64, bool)> = None;
    for &edge in trim_edges {
        let Some(edge_data) = brep.edges.get(edge) else {
            continue;
        };
        let he_a = edge_data.half_edge;
        let Some(he_a_data) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(he_b_data) = brep.half_edges.get(he_a_data.twin) else {
            continue;
        };
        // The support side of the trim: whichever half-edge's loop
        // face is not the strip.
        let f_a = brep.loops.get(he_a_data.loop_).map(|l| l.face);
        let f_b = brep.loops.get(he_b_data.loop_).map(|l| l.face);
        let support = match (f_a, f_b) {
            (Some(f), _) if f != strip_face => f,
            (_, Some(f)) if f != strip_face => f,
            _ => continue,
        };
        let Some(support_sof) = brep.faces.get(support).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let p0 = brep.vertices.get(he_a_data.origin).map(|v| v.position);
        let p1 = brep.vertices.get(he_b_data.origin).map(|v| v.position);
        let (Some(p0), Some(p1)) = (p0, p1) else {
            continue;
        };
        let sample = popkern_core::geometry::Point3::new(
            (p0.x + p1.x) * 0.5,
            (p0.y + p1.y) * 0.5,
            (p0.z + p1.z) * 0.5,
        );
        let normal_at = |geometry: &popkern_core::geometry::surface::SurfaceGeometry,
                         p: &popkern_core::geometry::Point3| {
            let eval = geometry as &dyn SurfaceEval;
            let (u, v) = eval.closest_point(p)?;
            eval.normal_at(u, v).normalized()
        };
        let Some(strip_normal) = normal_at(&strip_geometry, &sample) else {
            continue;
        };
        let Some(support_natural) = normal_at(&support_sof.geometry, &sample) else {
            continue;
        };
        let support_flag = if support_sof.orientation { 1.0 } else { -1.0 };
        let witness = strip_normal.dot(&(support_natural * support_flag));
        if best.is_none_or(|(w, _)| witness.abs() > w) {
            best = Some((witness.abs(), witness > 0.0));
        }
    }

    let (magnitude, agree) = best?;
    if magnitude < TANGENCY_WITNESS_FLOOR {
        return None;
    }
    if let Some(sof) = brep
        .faces
        .get_mut(strip_face)
        .and_then(|f| f.surface.as_mut())
    {
        sof.orientation = agree;
    }
    Some(agree)
}
