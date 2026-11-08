//! How two bodies whose boundaries meet in no intersection curve lie:
//! one inside the other, or apart (PK-1124).
//!
//! With no records, the boundaries cross nowhere, so every point of one
//! body's boundary is on the same side of the other, except where the
//! two only touch (a tangent point, a tangent curve, a shared face).
//! Points at those touches classify `OnBoundary` and say nothing. A
//! body's own vertices and edge samples are the likeliest to sit on a
//! touch (the inscribed ball's poles and seam lie on the cube's faces,
//! where it touches them), so face-interior points come first, several
//! per face, and the vertices and loop samples after.
//!
//! The verdict needs agreement: the decisive probes (Inside or Outside)
//! must all say the same. Disagreement means the boundaries do cross
//! and the intersection stage missed it; no decisive probe at all
//! means the relation is unknown. Both decline with
//! [`BoolError::RelationUndecided`]: reading either as "apart" served
//! an empty intersection of a cube and its inscribed ball.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

use crate::classify::point_in_solid_final;
use crate::error::{BoolError, PointClassification};
use crate::postprocess::zero_thickness::on_body_boundary;

/// Decisive probes that settle the relation once they agree.
const AGREEING_DECISIVE: usize = 3;
/// Face-interior points taken per face.
const INTERIOR_PER_FACE: usize = 3;
/// The most probes classified before the verdict is taken as it stands.
const PROBE_BUDGET: usize = 48;

/// What the probes of one body against another found.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct ProbeTally {
    /// Probes classified.
    pub probes: usize,
    /// Probes inside the other body.
    pub inside: usize,
    /// Probes outside it.
    pub outside: usize,
    /// Probes on its boundary.
    pub on_boundary: usize,
}

impl ProbeTally {
    /// The decided relation: Inside or Outside when the decisive probes
    /// agree, `OnBoundary` when no probe was decisive, `None` when they
    /// disagree.
    pub(crate) fn verdict(&self) -> Option<PointClassification> {
        match (self.inside > 0, self.outside > 0) {
            (true, false) => Some(PointClassification::Inside),
            (false, true) => Some(PointClassification::Outside),
            (false, false) => Some(PointClassification::OnBoundary),
            (true, true) => None,
        }
    }

    /// The decline for an undecided relation of `operand`'s probes.
    pub(crate) fn undecided(&self, operand: usize) -> BoolError {
        BoolError::RelationUndecided {
            operand,
            probes: self.probes,
            inside: self.inside,
            outside: self.outside,
            on_boundary: self.on_boundary,
        }
    }
}

/// Classify probe points of `inner` against `outer` until
/// [`AGREEING_DECISIVE`] decisive probes agree, one disagrees, or the
/// probes (at most [`PROBE_BUDGET`]) run out.
pub(crate) fn probe_relation(
    inner: &BRep,
    outer: &BRep,
    outer_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Result<ProbeTally, BoolError> {
    let probes = probe_points(inner);
    if probes.is_empty() {
        return Err(BoolError::DegenerateInput(
            "solid has no vertices".to_string(),
        ));
    }
    let mut tally = ProbeTally::default();
    for p in probes.iter().take(PROBE_BUDGET) {
        // Boolean operands are FINAL bodies: a ray hit inside one of a
        // face's holes is a void and must not count (PK-722).
        tally.probes += 1;
        // The ray vote has no on-boundary verdict of its own (a point
        // on the other body's face votes either way), so a probe at a
        // touch is recognised geometrically first.
        if on_body_boundary(p, outer, outer_aabbs, tolerance) {
            tally.on_boundary += 1;
            continue;
        }
        match point_in_solid_final(p, outer, outer_aabbs, tolerance)? {
            PointClassification::Inside => tally.inside += 1,
            PointClassification::Outside => tally.outside += 1,
            PointClassification::OnBoundary => tally.on_boundary += 1,
        }
        if tally.inside > 0 && tally.outside > 0 {
            break;
        }
        if tally.inside.max(tally.outside) >= AGREEING_DECISIVE {
            break;
        }
    }
    Ok(tally)
}

/// Probe points of `body`: face-interior points first, round-robin over
/// the faces so a budget covers every face, then the vertices, then
/// outer-loop samples.
fn probe_points(body: &BRep) -> Vec<Point3> {
    let mut faces: Vec<FaceHandle> = body.faces.keys().collect();
    faces.sort();
    let per_face: Vec<Vec<Point3>> = faces
        .iter()
        .map(|&f| face_interior_points(body, f, INTERIOR_PER_FACE))
        .collect();
    let mut out = Vec::new();
    for k in 0..INTERIOR_PER_FACE {
        out.extend(per_face.iter().filter_map(|pts| pts.get(k).copied()));
    }
    let mut vertices: Vec<_> = body.vertices.iter().collect();
    vertices.sort_by_key(|(v, _)| *v);
    out.extend(vertices.iter().map(|(_, v)| v.position));
    for &f in &faces {
        out.extend(body.face_outer_loop_sample_positions(f, 3));
    }
    out
}

/// Up to `max` points strictly inside `face`'s trim: its outer loop's
/// uv centroid, then an 8 × 8 grid over the loop's uv bounds, in an
/// order that spreads them (the grid read from its middle outward).
pub(crate) fn face_interior_points(body: &BRep, face: FaceHandle, max: usize) -> Vec<Point3> {
    let Some(geometry) = body
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    else {
        return Vec::new();
    };
    let uvs: Vec<(f64, f64)> = body
        .face_outer_loop_sample_positions(face, 4)
        .iter()
        .filter_map(|p| project_point_to_uv(geometry, p))
        .collect();
    if uvs.is_empty() {
        return Vec::new();
    }
    let (mut lo, mut hi) = (
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for &(u, v) in &uvs {
        lo = (lo.0.min(u), lo.1.min(v));
        hi = (hi.0.max(u), hi.1.max(v));
    }
    // A face bounded by its seam alone (a whole sphere: one loop down
    // the seam and back) projects its loop to a line, a box of no
    // width: there the surface's own range bounds the face.
    let ((u0, u1), (v0, v1)) = geometry.parameter_range();
    if hi.0 - lo.0 <= f64::EPSILON * (1.0 + lo.0.abs()) && u0.is_finite() && u1.is_finite() {
        (lo.0, hi.0) = (u0, u1);
    }
    if hi.1 - lo.1 <= f64::EPSILON * (1.0 + lo.1.abs()) && v0.is_finite() && v1.is_finite() {
        (lo.1, hi.1) = (v0, v1);
    }
    let n = uvs.len() as f64;
    let centre = (
        uvs.iter().map(|p| p.0).sum::<f64>() / n,
        uvs.iter().map(|p| p.1).sum::<f64>() / n,
    );
    // Grid stations 1..8 of 8 in each direction, nearest the middle
    // first, then alternating outward: 4, 3, 5, 2, 6, 1, 7.
    const ORDER: [u32; 7] = [4, 3, 5, 2, 6, 1, 7];
    let grid = ORDER.iter().flat_map(|&i| {
        ORDER.iter().map(move |&j| {
            (
                lo.0 + (hi.0 - lo.0) * f64::from(i) / 8.0,
                lo.1 + (hi.1 - lo.1) * f64::from(j) / 8.0,
            )
        })
    });
    let mut out: Vec<Point3> = Vec::new();
    for (u, v) in std::iter::once(centre).chain(grid) {
        if out.len() >= max {
            break;
        }
        if body.face_uv_trim_containment(face, geometry, u, v) != Some(FaceUvContainment::Interior)
        {
            continue;
        }
        let p = geometry.point_at(u, v);
        if out.iter().all(|q| q.distance_to(&p) > 0.0) {
            out.push(p);
        }
    }
    out
}

#[cfg(test)]
#[path = "operation_relation_probes_tests.rs"]
mod tests;
