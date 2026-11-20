//! Contact candidates (W3b, tangency design §4.8): two operand faces
//! touching at an isolated point where they meet in no curve enter the
//! table as candidates on both faces. Identity is the table's: a contact
//! joins the nearest compatible vertex, lies on a record's or an edge's
//! interior within the reach of it (a pave there, `incidence`), and a
//! contact on no record becomes a vertex of the result on each face
//! holding it (`VertexTable::isolated_contacts`).

use popkern_core::topology::BRep;

use super::candidates::{
    CandidateVertex, Incidence, SourceKind, handle_bits, own_tolerance, precise_radius,
};
use crate::parallel::ff_parallel::ff_parallel_contacts::FaceContact;

/// The candidate of each contact: a precise source (the touch is found
/// by Newton on the gap's gradient, design R5.8), widened by its faces'
/// own tolerances and half its gap (the point lies on A's face, B's the
/// gap away).
pub(super) fn contact_candidates(
    (a, b): (&BRep, &BRep),
    contacts: &[FaceContact],
    linear: f64,
) -> Vec<CandidateVertex> {
    contacts
        .iter()
        .map(|c| {
            let face_tol =
                |brep: &BRep, f| own_tolerance(brep.faces.get(f).and_then(|d| d.tolerance), linear);
            let radius = (face_tol(a, c.face_a) + face_tol(b, c.face_b))
                .max(precise_radius(linear))
                + 0.5 * c.gap.abs();
            let p = c.point;
            CandidateVertex {
                position: p,
                radius,
                zone: radius,
                source: SourceKind::Contact,
                incidence: vec![
                    Incidence::Face {
                        operand: 0,
                        face: c.face_a,
                        uv: Some(c.uv_a),
                    },
                    Incidence::Face {
                        operand: 1,
                        face: c.face_b,
                        uv: Some(c.uv_b),
                    },
                ],
                key: (
                    SourceKind::Contact,
                    0,
                    handle_bits(c.face_a),
                    handle_bits(c.face_b),
                    p.x.to_bits() ^ p.y.to_bits().rotate_left(21) ^ p.z.to_bits().rotate_left(42),
                ),
            }
        })
        .collect()
}

/// Where operands meeting in no curve touch (tangency design §4.8, W3d):
/// the vertex table over their contacts and interfering edge-face hits
/// (every hit is an edge resting on a face, with no record to end),
/// each resolved group's anchor once. Empty when the table declines.
pub(crate) fn touches_without_records(
    (a, b): (&BRep, &BRep),
    ix: &crate::intersect::IntersectionResult,
    linear: f64,
) -> Vec<popkern_core::geometry::Point3> {
    if ix.contacts.is_empty() && ix.edge_face_hits.is_empty() {
        return Vec::new();
    }
    let hits = super::candidates::interfering_hits((a, b), ix, &[], linear);
    super::table::build((a, b), (&hits, &ix.contacts), &[], linear)
        .isolated_contacts((a, b), linear)
}

/// Whether the edge-face hit `hit` is its edge resting on the face, not
/// crossing it: the edge's signed distance from the face's surface has
/// one sign either side of the hit, read at steps growing from a
/// millionth of the edge's span until both sides lie past the linear
/// reach (inside it the side is noise). A crossing changes sides. A hit
/// at the edge's end (its vertex) is judged on its one side.
pub(crate) fn hit_rests_on_face(
    (a, b): (&BRep, &BRep),
    hit: &CandidateVertex,
    linear: f64,
) -> bool {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let operands = [a, b];
    let (Some(&Incidence::Edge { operand, edge, t }), Some(&Incidence::Face { face, uv, .. })) = (
        hit.incidence
            .iter()
            .find(|i| matches!(i, Incidence::Edge { .. })),
        hit.incidence
            .iter()
            .find(|i| matches!(i, Incidence::Face { .. })),
    ) else {
        return false;
    };
    let Some(curve) = operands[operand]
        .edges
        .get(edge)
        .and_then(|e| e.curve.as_ref())
    else {
        return false;
    };
    let Some(surface) = operands[1 - operand]
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    else {
        return false;
    };
    let reach = popkern_core::tolerance::Tolerance::default()
        .with_linear(linear)
        .linear_reach();
    let (lo, hi) = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    let near = uv.unwrap_or_else(|| {
        popkern_core::geometry::uv_projection::project_point_to_uv(surface, &hit.position)
            .unwrap_or((0.0, 0.0))
    });
    let signed = |s: f64| -> Option<f64> {
        let p = curve.geometry.point_at(s);
        let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv_near(
            surface, &p, near.0, near.1,
        )?;
        let n = surface.normal_at(u, v).normalized()?;
        Some((p - surface.point_at(u, v)).dot(&n))
    };
    let mut step = 1e-6 * (hi - lo);
    while step <= hi - lo {
        let (sm, sp) = (t - step, t + step);
        let side = |s: f64| (lo..=hi).contains(&s).then(|| signed(s)).flatten();
        match (side(sm), side(sp)) {
            (Some(dm), Some(dp)) if dm.abs() > reach && dp.abs() > reach => {
                return dm.signum() == dp.signum();
            }
            // At the edge's end, its vertex: the edge leaves the face on
            // one side, and an edge of the vertex crossing it would end a
            // record there; on no record the vertex rests on the face.
            (None, Some(d)) | (Some(d), None) if d.abs() > reach => return true,
            (None, None) => return false,
            _ => {}
        }
        step *= 4.0;
    }
    false
}

#[cfg(test)]
#[path = "contact_candidates_tests.rs"]
mod tests;
