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

use super::ids::FaceId;
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

/// A point on `face` a short way in from its first fin: the fin's
/// midpoint stepped along the inward direction (normal × fin tangent)
/// by 1% of the fin, at least 100 tolerances, and projected back onto
/// the surface.
pub(crate) fn inner_point(
    c: &CellComplex,
    operands: [&BRep; 2],
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
    let fin = c.loops.get(f.loops.first()?.index())?.first;
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
    let on_face = mid + into * (len * 1e-2).max(linear * 100.0);
    let (u2, v2) = surf.geometry.closest_point(&on_face)?;
    Some(InnerPoint {
        point: surf.geometry.point_at(u2, v2),
        uv: (u2, v2),
        fin_length: len,
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
    operands: [&BRep; 2],
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
