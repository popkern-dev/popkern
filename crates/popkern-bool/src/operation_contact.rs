//! Bodies that only touch (PK-756).
//!
//! Two solids meeting on coincident faces whose normals oppose, and
//! otherwise at most at points, have disjoint interiors: near each
//! contact face the bodies lie on opposite sides of it, no surface of
//! one crosses the other (every record is a point or runs along a
//! contact face's rim, where the contact ends), and containment would
//! need the contact faces cooriented.
//! So their intersection is empty and A minus B is A.
//!
//! The full pipeline got this wrong on the cycloid module's output pin
//! standing on its flange: the pin's top circle met the flange face as
//! a zero-length record, splitting a zero-height sliver off the pin
//! wall, sampled ON the flange plane and voted Inside; the cellular
//! assembly declined (a valence-1 vertex) and the graph-merge fallback
//! served the two coincident end discs as a 301 mm³ "overlap" that the
//! interference check reported as a collision.
//!
//! Coincidence detection pairs faces by SURFACE, not by region, so two
//! cooriented faces that only abut on a shared plane also arrive as an
//! overlap: the cycloid M module's flange hub (bottom annulus r 27.03 to
//! 29.03) beside its bearing's inner race (bottom annulus r 29.03 to
//! 33.9) on one plane. Such a pair shares no area, so it says nothing
//! about the interiors and is set aside; the opposite contact (hub wall
//! on race bore) still decides. Without that, the full pipeline built
//! the contact cylinder as a two-face "overlap" of 16587 mm³.
//!
//! The union of touching bodies merges their contact faces, so it keeps
//! the full pipeline.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::topology::FaceHandle;
use popkern_core::validation::FaceUvContainment;

use crate::coincidence::types::CoincidenceOrientation;
use crate::error::BoolOperation;
use crate::intersect::IntersectionResult;

/// Samples per record when measuring its length.
const LENGTH_SAMPLES: usize = 8;

/// The result when `a` and `b` only touch (see the module docs):
/// `Some` for an intersection (empty) or a subtraction (`a`), `None`
/// otherwise — a record crossing anywhere but a contact face's rim, a
/// cooriented coincidence that shares area, no opposite contact at
/// all, or a union.
pub(crate) fn touching_only(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    intersection: &IntersectionResult,
    tolerance: &Tolerance,
) -> Option<BRep> {
    if matches!(op, BoolOperation::Union) {
        return None;
    }
    if intersection.coincidence.face_overlaps.is_empty() {
        return None;
    }
    // Cooriented pairs that only abut are no contact (module docs); any
    // other cooriented pair shares area, and the bodies may overlap.
    let mut overlaps = Vec::new();
    for o in &intersection.coincidence.face_overlaps {
        if o.orientation == CoincidenceOrientation::Opposite {
            overlaps.push(o);
        } else if !abutting(a, o.face_a, b, o.face_b) {
            trace_decline(&format!(
                "cooriented A:{:?} B:{:?} shares area",
                o.face_a, o.face_b
            ));
            return None;
        }
    }
    // An opposite contact is what rules containment out.
    if overlaps.is_empty() {
        trace_decline("no opposite contact");
        return None;
    }
    // Every record must be a point or run along the rim of a contact
    // face (the pin wall meets the flange plane only along its own top
    // edge): such a record is where the contact ends, not a crossing.
    let rims: Vec<&popkern_core::geometry::curve::CurveOnEdge> = overlaps
        .iter()
        .flat_map(|o| [(a, o.face_a), (b, o.face_b)])
        .flat_map(|(brep, face)| face_edges(brep, face))
        .collect();
    let near = tolerance.linear() * 10.0;
    let on_rim = |p: &popkern_core::geometry::Point3| {
        rims.iter().any(|c| {
            let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
            c.geometry
                .closest_point(p)
                .map(|t| t.clamp(lo, hi))
                .is_some_and(|t| c.geometry.point_at(t).distance_to(p) <= near)
        })
    };
    // A trim residue (`residue_owner_a`) is the other operand's surface
    // crossing ONE face beyond the other face's trim: no boundary of
    // the second body is there, so it is no crossing of the bodies
    // (PK-847: the six-axis arm's wrist casting's ring and web planes
    // left residues along a butted tube's whole wall, and the contact
    // classified as a crossing).
    let contact_only = intersection.records.iter().all(|r| {
        if r.residue_owner_a.is_some() {
            return true;
        }
        let c = &r.curve.curve_3d;
        let (t0, t1) = c.parameter_range();
        let samples: Vec<_> = (0..=LENGTH_SAMPLES)
            .map(|i| c.point_at(t0 + (t1 - t0) * i as f64 / LENGTH_SAMPLES as f64))
            .collect();
        let length: f64 = samples.windows(2).map(|w| w[0].distance_to(&w[1])).sum();
        length <= near || samples.iter().all(on_rim)
    });
    if !contact_only {
        trace_decline("a record crosses off the contact rims");
        return None;
    }
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [exec] -> touching only ({} opposite contact face(s), {} point record(s))",
            overlaps.len(),
            intersection.records.len()
        );
    }
    match op {
        BoolOperation::Intersection => Some(BRep::new(*tolerance)),
        BoolOperation::Subtraction => Some(a.clone()),
        BoolOperation::Union => None,
    }
}

fn trace_decline(why: &str) {
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!("  [exec] touching-only declined: {why}");
    }
}

/// Whether coincident faces `fa` of `a` and `fb` of `b` share no area:
/// neither has an interior point inside the other's trim. Exact when
/// their boundaries do not cross, which `touching_only` establishes
/// before the verdict is used (every record lies on a contact rim);
/// then the regions are disjoint, nested, or equal, and a nested or
/// equal pair puts an interior point of one inside the other. No
/// interior point or no projection counts as sharing (decline).
fn abutting(a: &BRep, fa: FaceHandle, b: &BRep, fb: FaceHandle) -> bool {
    !reaches_into(a, fa, b, fb) && !reaches_into(b, fb, a, fa)
}

/// Whether an interior point of `fx` lies on `fy`'s region (inside or
/// on its trim boundary), or no verdict could be reached.
fn reaches_into(x: &BRep, fx: FaceHandle, y: &BRep, fy: FaceHandle) -> bool {
    let Some(p) = face_interior_point(x, fx) else {
        return true;
    };
    let Some(surface) = y.faces.get(fy).and_then(|f| f.surface.as_ref()) else {
        return true;
    };
    let Some((u, v)) = project_point_to_uv(&surface.geometry, &p) else {
        return true;
    };
    y.face_uv_trim_containment(fy, &surface.geometry, u, v) != Some(FaceUvContainment::Outside)
}

/// Steps inward from the outer loop, as fractions of the face's uv
/// extent: thin faces (a 2 mm annulus on a 58 mm disc) need the small
/// ones, which the lineage sampler's centroid-and-grid search misses.
const INWARD_STEPS: [f64; 3] = [1e-2, 1e-3, 1e-4];

/// A point strictly inside `face`'s trim: the lineage sampler's pick,
/// else a point just inside the outer loop, beside the midpoint of a
/// sampled boundary chord.
fn face_interior_point(brep: &BRep, face: FaceHandle) -> Option<popkern_core::geometry::Point3> {
    if let Some((p, _)) = crate::entity_lineage::interior_point(brep, face) {
        return Some(p);
    }
    let geometry = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let uvs: Vec<(f64, f64)> = brep
        .face_outer_loop_sample_positions(face, 8)
        .iter()
        .filter_map(|p| project_point_to_uv(geometry, p))
        .collect();
    let (mut lo, mut hi) = (
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for &(u, v) in &uvs {
        lo = (lo.0.min(u), lo.1.min(v));
        hi = (hi.0.max(u), hi.1.max(v));
    }
    let extent = (hi.0 - lo.0).max(hi.1 - lo.1);
    if extent.is_nan() || extent <= 0.0 {
        return None;
    }
    for step in INWARD_STEPS {
        for w in uvs.windows(2) {
            let ((u0, v0), (u1, v1)) = (w[0], w[1]);
            let (du, dv) = (u1 - u0, v1 - v0);
            let len = (du * du + dv * dv).sqrt();
            if len.is_nan() || len <= 0.0 {
                continue;
            }
            let (mu, mv) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
            // The interior is on one side of the chord; try both.
            for side in [1.0, -1.0] {
                let (u, v) = (
                    mu - side * dv / len * step * extent,
                    mv + side * du / len * step * extent,
                );
                if brep.face_uv_trim_containment(face, geometry, u, v)
                    == Some(FaceUvContainment::Interior)
                {
                    return Some(geometry.point_at(u, v));
                }
            }
        }
    }
    None
}

/// The edge curves bounding `face` (every loop).
fn face_edges(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> Vec<&popkern_core::geometry::curve::CurveOnEdge> {
    let mut out = Vec::new();
    let Some(f) = brep.faces.get(face) else {
        return out;
    };
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(c) = brep.edges.get(h.edge).and_then(|e| e.curve.as_ref()) {
                out.push(c);
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    out
}
