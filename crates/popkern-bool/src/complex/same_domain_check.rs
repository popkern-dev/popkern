//! Same-domain completeness (plan §4.4): every fragment lying in a
//! coincident patch is grouped with its partner's copy.
//!
//! The walk cuts a coincident patch once on each operand's face, and
//! `same_domain_groups` joins the two copies into one complex face. A
//! copy left alone (the walks split the patch differently, or one walk
//! dropped its piece) lifts as a face carried by one operand that lies
//! on its partner face: two coincident faces where the complex should
//! hold one, or one where it should hold two owners. Radial ordering
//! would then see a tie no ladder can break, or membership would read
//! a skin as a wall. `SameDomainMismatch` reports it at the join,
//! naming the face, its owner and the partner it lies on.
//!
//! The test is geometric and per face: a point just inside the face
//! (in from its first fin, as membership probes it) lies on the partner
//! face's surface within tolerance and inside its trim region. A
//! region test with no verdict (a periodic chart) reports nothing.

use std::collections::BTreeSet;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

use super::ids::{FaceId, FinId};
use super::model::{CellComplex, FaceOwner};

/// A face carried by one operand that lies on its coincident partner.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SameDomainMismatch {
    /// The complex face.
    pub face: FaceId,
    /// Its only owner.
    pub owner: FaceOwner,
    /// The partner face (of the other operand) it lies on.
    pub partner: FaceHandle,
    /// The fragments the face lifts.
    pub fragments: Vec<usize>,
}

/// A point just inside a complex face, on its primary owner's surface.
#[derive(Debug, Clone, Copy)]
pub(crate) struct InnerPoint {
    /// The point.
    pub point: Point3,
    /// Its parameters on the owner's surface.
    pub uv: (f64, f64),
    /// The length of the fin it was stepped in from.
    pub fin_length: f64,
}

/// A point on `face` a short way in from one of its fins: the fin's
/// midpoint stepped along the inward direction (normal × fin tangent)
/// and projected back onto the surface. The step is 1% of the fin, at
/// least 100 tolerances, and a candidate stands only when every other
/// edge of the face is farther from it than the step: a step that left
/// the face crossed another edge, which then lies nearer than the step
/// (measured: a stadium's rod-cap horn 7.5e-4 thick, stepped 9.5e-4 in
/// from its rim arc, landed on the box cap beyond it). A rejected
/// candidate shrinks the step, then tries the face's other fins; `None`
/// when none stands.
pub(crate) fn inner_point(
    c: &CellComplex,
    operands: &[&BRep],
    linear: f64,
    face: FaceId,
) -> Option<InnerPoint> {
    let f = c.faces.get(face.index())?;
    let owner = f.owners.first()?;
    let surf = operands[owner.operand]
        .faces
        .get(owner.face)?
        .surface
        .as_ref()?;
    let fins: Vec<FinId> = f
        .loops
        .iter()
        .flat_map(|&l| c.loop_fins(l))
        .take(16)
        .collect();
    for shrink in [1.0, 0.25, 1.0 / 16.0] {
        for &fin in &fins {
            let Some(at) = step_in(c, surf, linear, fin, shrink) else {
                continue;
            };
            if clear_of_other_edges(c, &f.loops, fin, &at) {
                return Some(at.inner);
            }
        }
    }
    None
}

/// A candidate inner point and the step that reached it.
struct Stepped {
    inner: InnerPoint,
    step: f64,
}

fn step_in(
    c: &CellComplex,
    surf: &popkern_core::geometry::SurfaceOnFace,
    linear: f64,
    fin: FinId,
    shrink: f64,
) -> Option<Stepped> {
    let edge = &c.edges[c.fins[fin.index()].edge.index()];
    let (pa, pb) = (
        c.vertices[edge.start.index()].position,
        c.vertices[edge.end.index()].position,
    );
    let mid = match &edge.curve {
        Some(cv) => {
            use popkern_core::geometry::curve::CurveEval;
            cv.geometry.point_at(0.5 * (cv.t_start + cv.t_end))
        }
        None => pa + (pb - pa) * 0.5,
    };
    let (u, v) = surf.geometry.closest_point(&mid)?;
    let n = surf.folded_normal(u, v)?;
    let t = (pb - pa).normalized().or_else(|| {
        edge.curve.as_ref().and_then(|cv| {
            use popkern_core::geometry::curve::CurveEval;
            cv.geometry
                .derivative_at(0.5 * (cv.t_start + cv.t_end), 1)
                .normalized()
        })
    })?;
    let tf = if c.fins[fin.index()].forward {
        t
    } else {
        t * -1.0
    };
    let len = pa.distance_to(&pb).max(linear * 1e3);
    let into = n.cross(&tf).normalized()?;
    let step = (len * 1e-2).max(linear * 100.0) * shrink;
    if step < linear * 10.0 {
        return None;
    }
    let on_face = mid + into * step;
    let (u2, v2) = surf.geometry.closest_point(&on_face)?;
    Some(Stepped {
        inner: InnerPoint {
            point: surf.geometry.point_at(u2, v2),
            uv: (u2, v2),
            fin_length: len,
        },
        step,
    })
}

/// Whether every edge of `loops` other than `fin`'s lies farther than
/// the step from the candidate.
fn clear_of_other_edges(
    c: &CellComplex,
    loops: &[super::ids::LoopId],
    fin: FinId,
    at: &Stepped,
) -> bool {
    let own = c.fins[fin.index()].edge;
    let p = at.inner.point;
    loops.iter().flat_map(|&l| c.loop_fins(l)).all(|g| {
        let e = c.fins[g.index()].edge;
        if e == own {
            return true;
        }
        let edge = &c.edges[e.index()];
        let (a, b) = (
            c.vertices[edge.start.index()].position,
            c.vertices[edge.end.index()].position,
        );
        let d = match &edge.curve {
            Some(cv) => {
                use popkern_core::geometry::curve::CurveEval;
                let t = popkern_core::geometry::curve::closest_point_in_window(
                    &cv.geometry,
                    &p,
                    cv.t_start,
                    cv.t_end,
                );
                cv.geometry.point_at(t).distance_to(&p)
            }
            None => {
                let ab = b - a;
                let k = ((p - a).dot(&ab) / ab.dot(&ab).max(1e-300)).clamp(0.0, 1.0);
                (a + ab * k).distance_to(&p)
            }
        };
        d > at.step
    })
}

/// Whether `p` lies on `face` of `brep`: within `reach` of its surface
/// and inside its trim region. `None` when the region test gives no
/// verdict.
fn lies_on_face(brep: &BRep, face: FaceHandle, p: &Point3, reach: f64) -> Option<bool> {
    let surf = brep.faces.get(face)?.surface.as_ref()?;
    let (u, v) = surf.geometry.closest_point(p)?;
    if surf.geometry.point_at(u, v).distance_to(p) > reach {
        return Some(false);
    }
    crate::classify_uv_loops::point_in_face_region_exact(brep, face, u, v)
}

/// Every face of `c` carried by one operand that lies on that owner's
/// coincident partner (`overlaps`: (A face, B face) pairs).
pub(crate) fn same_domain_mismatches(
    c: &CellComplex,
    operands: &[&BRep],
    overlaps: &[(FaceHandle, FaceHandle)],
    linear: f64,
) -> Vec<SameDomainMismatch> {
    if overlaps.is_empty() {
        return Vec::new();
    }
    let pairs: BTreeSet<(FaceHandle, FaceHandle)> = overlaps.iter().copied().collect();
    let mut out = Vec::new();
    for (k, f) in c.faces.iter().enumerate() {
        let [owner] = f.owners.as_slice() else {
            continue;
        };
        let partners = pairs.iter().filter_map(|&(a, b)| match owner.operand {
            0 if a == owner.face => Some(b),
            1 if b == owner.face => Some(a),
            _ => None,
        });
        let face = FaceId::from_index(k);
        let mut probe = None;
        for partner in partners {
            let Some(at) = probe.get_or_insert_with(|| inner_point(c, operands, linear, face))
            else {
                break;
            };
            if lies_on_face(
                operands[1 - owner.operand],
                partner,
                &at.point,
                linear * 10.0,
            ) == Some(true)
            {
                out.push(SameDomainMismatch {
                    face,
                    owner: *owner,
                    partner,
                    fragments: f.source_fragments.clone(),
                });
                break;
            }
        }
    }
    out
}

#[cfg(test)]
#[path = "same_domain_check_tests.rs"]
mod tests;
