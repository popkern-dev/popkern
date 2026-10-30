//! Shallow-crossing junction pins.
//!
//! A record leaving a face's trim at a shallow angle θ to the boundary
//! edge has its uv crossing displaced ALONG the record by the uv fit's
//! lateral slop over sin θ. The exact junction — the boundary edge's
//! pierce of the other operand's surface — is in the pierce inventory,
//! but the exact-match pin (`clip_windows`) matches window bounds to
//! junctions only within span × 1e-6 and only once the pierce windows
//! cleared the boundary. A piece end off its junction by more meets
//! nothing there: the other face's walk sees a dangling end, and the
//! record is used by one face only (measured: a boss cap 1e-3 past a
//! ring's extreme cuts the web's bottom face along a line leaving the
//! ring's rim arc at 0.34°; the clipped end lay 1.6e-4 along the line
//! from the junction, and the cap never split along the piece).
//!
//! A piece end pins onto a junction when their distance is what the
//! crossing angle explains: within 10 × linear / sin θ, θ measured
//! between the record and the boundary edge through the junction, and
//! within the exact-match pin's own 3D gate. A steep crossing explains
//! almost nothing, so only shallow crossings move.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_pierce_windows::{PierceInventoryCache, pierce_crossing_params};

/// The unit tangent of the boundary edge of `face` passing through
/// `p` (within 10 × linear), if any.
fn boundary_tangent_at(brep: &BRep, face: FaceHandle, p: &Point3) -> Option<[f64; 3]> {
    let band = 10.0 * brep.tolerance.linear();
    let mut best: Option<(f64, [f64; 3])> = None;
    for eh in crate::clip_endpoint_refine::trim_boundary_edges(brep, face) {
        let Some(c) = brep.edges.get(eh).and_then(|e| e.curve.as_ref()) else {
            continue;
        };
        let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
        let t = popkern_core::geometry::curve::closest_point_in_window(&c.geometry, p, lo, hi);
        let d = c.geometry.point_at(t).distance_to(p);
        if d >= band || best.as_ref().is_some_and(|(bd, _)| *bd <= d) {
            continue;
        }
        let v = c.geometry.tangent_at(t);
        let n = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
        if n > 0.0 {
            best = Some((d, [v.x / n, v.y / n, v.z / n]));
        }
    }
    best.map(|(_, v)| v)
}

/// Pin `piece`'s ends onto exact junctions their offset is explained
/// by the crossing angle (see the module docs). A face without a
/// surface leaves the piece as it is.
pub(super) fn pin_shallow_crossing_ends(
    piece: &mut CurveOnSurface,
    (brep_a, face_a): (&BRep, FaceHandle),
    (brep_b, face_b): (&BRep, FaceHandle),
    record: &CurveOnSurface,
    surfaces: (Option<&SurfaceGeometry>, Option<&SurfaceGeometry>),
    cache: &mut PierceInventoryCache,
) {
    let (Some(surface_a), Some(surface_b)) = surfaces else {
        return;
    };
    let linear = brep_a.tolerance.linear().max(brep_b.tolerance.linear());
    let gate = 2000.0 * linear;
    let mut junctions: Vec<(Point3, FaceHandle, bool)> = Vec::new();
    for (bx, fx, ob, of, x_is_a) in [
        (brep_a, face_a, brep_b, face_b, true),
        (brep_b, face_b, brep_a, face_a, false),
    ] {
        for (_, p) in
            pierce_crossing_params(bx, fx, ob, of, record, cache, x_is_a).unwrap_or_default()
        {
            junctions.push((p, fx, x_is_a));
        }
    }
    if junctions.is_empty() {
        return;
    }
    let (t0, t1) = piece.curve_3d.parameter_range();
    for (end, t) in [(false, t0), (true, t1)] {
        let p_end = piece.curve_3d.point_at(t);
        let tr = piece.curve_3d.tangent_at(t);
        let tr_len = (tr.x * tr.x + tr.y * tr.y + tr.z * tr.z).sqrt();
        if tr_len <= 0.0 {
            continue;
        }
        let mut best: Option<(f64, Point3)> = None;
        for &(pj, fx, x_is_a) in &junctions {
            let d = p_end.distance_to(&pj);
            if d >= gate || d <= linear || best.as_ref().is_some_and(|(bd, _)| *bd <= d) {
                continue;
            }
            let brep = if x_is_a { brep_a } else { brep_b };
            let Some(te) = boundary_tangent_at(brep, fx, &pj) else {
                continue;
            };
            let (rx, ry, rz) = (tr.x / tr_len, tr.y / tr_len, tr.z / tr_len);
            let cross = (
                ry * te[2] - rz * te[1],
                rz * te[0] - rx * te[2],
                rx * te[1] - ry * te[0],
            );
            let sin = (cross.0 * cross.0 + cross.1 * cross.1 + cross.2 * cross.2).sqrt();
            if sin <= 0.0 || d > 10.0 * linear / sin {
                continue;
            }
            best = Some((d, pj));
        }
        if let Some((_, pj)) = best {
            crate::clip_endpoint_refine::set_endpoint_3d_uv(
                piece,
                end,
                pj,
                surface_a,
                surface_b,
                &brep_a.tolerance,
            );
        }
    }
}
