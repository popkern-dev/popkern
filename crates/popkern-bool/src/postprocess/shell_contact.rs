//! Contact between different shells of a served result (PK-964, D1).
//!
//! Operands that touch without sharing any topology — a ball resting on
//! a face, a rod lying on a block, an edge standing on a face — come out
//! of a union as separate shells that merely touch. Every manifold
//! check passes such a body (each shell is closed and clean on its
//! own), yet it is non-manifold in fact: the shells meet with zero
//! thickness. The intersector does not always report the contact (a
//! tangency is a double root it may not see), so the test here is
//! geometric and needs no intersection evidence.
//!
//! For every pair of faces on different shells whose boxes come within
//! the linear tolerance, seeds sampled on both faces (vertices, loop
//! samples, an interior grid) are refined by alternating closest-point
//! projection between the two surfaces, which converges to a local
//! minimum of their distance. A minimum within the linear tolerance
//! whose feet both lie inside or on their faces' trims is a contact.
//! The face-pair box test prunes every pair that cannot touch, so a
//! cavity shell with wall thickness pays nothing.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

use crate::error::{ContactKind, NonManifoldContact};

/// Boundary samples per edge of a face's outer loop.
const LOOP_SAMPLES: usize = 6;
/// Interior grid resolution per uv direction.
const GRID: usize = 6;
/// Seeds refined per face pair (the closest by initial distance).
const REFINED_SEEDS: usize = 4;
/// Alternating-projection rounds per seed.
const ROUNDS: usize = 24;

/// The closest point of `face`'s surface to `p`, with its uv.
fn foot(surface: &SurfaceGeometry, p: &Point3) -> Option<(Point3, (f64, f64))> {
    let (u, v) = project_point_to_uv(surface, p)?;
    Some((surface.point_at(u, v), (u, v)))
}

fn surface_of(brep: &BRep, face: FaceHandle) -> Option<&SurfaceGeometry> {
    Some(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}

fn within_trim(brep: &BRep, face: FaceHandle, surface: &SurfaceGeometry, uv: (f64, f64)) -> bool {
    brep.face_uv_trim_containment(face, surface, uv.0, uv.1)
        .is_some_and(|c| c != FaceUvContainment::Outside)
}

/// Seed points on `face`: its loop samples and an interior uv grid.
fn seeds(brep: &BRep, face: FaceHandle, surface: &SurfaceGeometry) -> Vec<Point3> {
    let mut out = brep.face_outer_loop_sample_positions(face, LOOP_SAMPLES);
    let uvs: Vec<(f64, f64)> = out
        .iter()
        .filter_map(|p| project_point_to_uv(surface, p))
        .collect();
    let (mut lo, mut hi) = (
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for &(u, v) in &uvs {
        lo = (lo.0.min(u), lo.1.min(v));
        hi = (hi.0.max(u), hi.1.max(v));
    }
    if lo.0.is_finite() && hi.0.is_finite() {
        for i in 0..GRID {
            for j in 0..GRID {
                let u = lo.0 + (hi.0 - lo.0) * (i as f64 + 0.5) / GRID as f64;
                let v = lo.1 + (hi.1 - lo.1) * (j as f64 + 0.5) / GRID as f64;
                if brep.face_uv_trim_containment(face, surface, u, v)
                    == Some(FaceUvContainment::Interior)
                {
                    out.push(surface.point_at(u, v));
                }
            }
        }
    }
    if let Some((p, _)) = crate::entity_lineage::interior_point(brep, face) {
        out.push(p);
    }
    out
}

/// How far a contact extends: a point, a curve, or an area.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ContactExtent {
    /// The faces meet at an isolated point.
    Point,
    /// The faces meet along a curve through the point.
    Line,
    /// The faces coincide over a neighbourhood of the point.
    Area,
}

/// A touch between faces of two different shells.
#[derive(Debug, Clone)]
pub(crate) struct ShellFinding {
    /// Where, and what kind (a line contact reports as an edge).
    pub contact: NonManifoldContact,
    /// The touching faces.
    pub faces: (FaceHandle, FaceHandle),
    /// How far the touch extends.
    pub extent: ContactExtent,
}

/// Probe directions around a contact, every 45° in the tangent plane.
const EXTENT_DIRECTIONS: usize = 8;

/// The extent of the touch of `f` and `g` at `at`: how many directions
/// around it in `f`'s tangent plane stay on both faces. None is a
/// point; one opposite pair is a line through it; more is an area.
fn contact_extent(
    brep: &BRep,
    (f, sf): (FaceHandle, &SurfaceGeometry),
    (g, sg): (FaceHandle, &SurfaceGeometry),
    at: Point3,
    linear: f64,
) -> ContactExtent {
    let Some((u, v)) = project_point_to_uv(sf, &at) else {
        return ContactExtent::Point;
    };
    let normal = sf.normal_at(u, v);
    let (t1, t2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(normal);
    let size = crate::intersect::compute_face_aabbs(brep)
        .into_iter()
        .find(|(h, _)| *h == f)
        .map_or(1.0, |(_, b)| b.diagonal());
    // At the face's own scale: an area touch coincides over a
    // neighbourhood that size, while a tangency separates there by the
    // curvature term (a 5 mm rod 1e-4 off its contact line sits 1e-6
    // off the plane: one tolerance, not a separation).
    let step = (5.0e-2 * size).max(100.0 * linear);
    let on_both: Vec<bool> = (0..EXTENT_DIRECTIONS)
        .map(|k| {
            let angle = std::f64::consts::TAU * k as f64 / EXTENT_DIRECTIONS as f64;
            use popkern_core::geometry::PortableMath;
            let d = t1 * angle.pcos() + t2 * angle.psin();
            let probe = at + d * step;
            let Some((p, uv_f)) = foot(sf, &probe) else {
                return false;
            };
            let Some((q, uv_g)) = foot(sg, &p) else {
                return false;
            };
            p.distance_to(&q) <= linear
                && within_trim(brep, f, sf, uv_f)
                && within_trim(brep, g, sg, uv_g)
        })
        .collect();
    let count = on_both.iter().filter(|&&b| b).count();
    let half = EXTENT_DIRECTIONS / 2;
    let opposite_pair = (0..half).any(|k| on_both[k] && on_both[k + half]);
    match count {
        0 => ContactExtent::Point,
        1 | 2 if count == 1 || opposite_pair => ContactExtent::Line,
        _ => ContactExtent::Area,
    }
}

/// The contact point of faces `f` and `g`, if they touch.
fn face_contact(brep: &BRep, f: FaceHandle, g: FaceHandle, linear: f64) -> Option<Point3> {
    let (sf, sg) = (surface_of(brep, f)?, surface_of(brep, g)?);
    let mut starts: Vec<(f64, Point3, Point3)> = Vec::new();
    for p in seeds(brep, f, sf) {
        if let Some((q, _)) = foot(sg, &p) {
            starts.push((p.distance_to(&q), p, q));
        }
    }
    for q in seeds(brep, g, sg) {
        if let Some((p, _)) = foot(sf, &q) {
            starts.push((p.distance_to(&q), p, q));
        }
    }
    starts.sort_by(|a, b| a.0.total_cmp(&b.0));
    for &(_, p0, q0) in starts.iter().take(REFINED_SEEDS) {
        let (mut p, mut q) = (p0, q0);
        let (mut uv_f, mut uv_g) = (project_point_to_uv(sf, &p)?, project_point_to_uv(sg, &q)?);
        for _ in 0..ROUNDS {
            let (np, nuv_f) = foot(sf, &q)?;
            let (nq, nuv_g) = foot(sg, &np)?;
            let moved = np.distance_to(&p).max(nq.distance_to(&q));
            (p, q, uv_f, uv_g) = (np, nq, nuv_f, nuv_g);
            if moved < 0.1 * linear {
                break;
            }
        }
        if p.distance_to(&q) <= linear
            && within_trim(brep, f, sf, uv_f)
            && within_trim(brep, g, sg, uv_g)
        {
            return Some(p + (q - p) * 0.5);
        }
    }
    None
}

/// Every place two different shells of `brep` touch, in face-handle
/// order, one per location.
pub(crate) fn shell_contacts(brep: &BRep, tolerance: &Tolerance) -> Vec<NonManifoldContact> {
    let linear = tolerance.linear();
    let mut out: Vec<NonManifoldContact> = Vec::new();
    for finding in shell_findings(brep, tolerance) {
        let at = finding.contact.location;
        if out
            .iter()
            .all(|c| c.location.distance_to(&at) > 10.0 * linear)
        {
            out.push(finding.contact);
        }
    }
    out
}

/// Every touching face pair of different shells, with where and how far
/// the pair touches (one finding per pair: a contact record must name
/// every pair it excuses).
pub(crate) fn shell_findings(brep: &BRep, tolerance: &Tolerance) -> Vec<ShellFinding> {
    if brep.back_shell_count() < 2 {
        return Vec::new();
    }
    let linear = tolerance.linear();
    let mut faces: Vec<(FaceHandle, Aabb)> = crate::intersect::compute_face_aabbs(brep);
    faces.sort_by_key(|(f, _)| *f);
    let shell = |f: FaceHandle| brep.faces.get(f).map(|face| face.shell);
    let mut out: Vec<ShellFinding> = Vec::new();
    for (i, (f, box_f)) in faces.iter().enumerate() {
        for (g, box_g) in &faces[i + 1..] {
            if shell(*f) == shell(*g) || !box_f.overlaps(box_g, linear) {
                continue;
            }
            let Some(at) = face_contact(brep, *f, *g, linear) else {
                continue;
            };
            let (Some(sf), Some(sg)) = (surface_of(brep, *f), surface_of(brep, *g)) else {
                continue;
            };
            let extent = contact_extent(brep, (*f, sf), (*g, sg), at, linear);
            {
                out.push(ShellFinding {
                    contact: NonManifoldContact {
                        kind: if extent == ContactExtent::Point {
                            ContactKind::Vertex
                        } else {
                            ContactKind::Edge
                        },
                        location: at,
                    },
                    faces: (*f, *g),
                    extent,
                });
            }
        }
    }
    out
}

#[cfg(test)]
#[path = "shell_contact_tests.rs"]
mod tests;
