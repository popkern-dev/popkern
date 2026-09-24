//! Ring sense of a boundary-ring (pushed) vertex-blend boundary,
//! fixed by the solid rather than by an arc's parameter direction
//! (PK-517).
//!
//! `reorient_arcs_for_hex_closure` builds a closed cycle of fillet
//! end arcs, but the successor walk fixes the SENSE of that cycle
//! from arc 0's parameter direction, and the angular sort fixes it
//! from a best-fit axis whose sign is arbitrary. Either way the
//! cycle's rotational sense is an accident of how the extractor
//! parameterized one arc. The boundary-ring rung re-extracts every
//! end arc at pushed stations, so the same corner can arrive with
//! either sense — measured on the carved pentagonal pyramid: the
//! unpushed boundary circulated one way, the pushed ring the other,
//! and every cap builder on the ring assembled a mirrored cap
//! (cap-contains-curve 1.6e-1, G1 normal deficit 0.999) that the
//! acceptance gates then declined.
//!
//! The ring carve BUILDS the blend face from this ring, so the ring
//! must circulate in that face's loop sense. The sense is read
//! locally and exactly at every bridge rather than from an averaged
//! outward direction (a saddle's support normals sum to the concave
//! side and mislead any average): bridge `k` runs across support
//! face `F` from arc `k`'s end to arc `k+1`'s start, and the carve
//! hands the vertex side of that cut to the cap. With `n_F` the
//! face-sensed outward normal of `F`, `t` the bridge direction and
//! `c` the in-surface direction from the bridge toward the vertex,
//! `(n_F × t) · c` has one sign for every bridge of a ring in loop
//! sense and the other for every bridge of a mirrored ring. The
//! bridges vote; a unanimous vote decides. A split or empty vote
//! (no support face found, degenerate bridge) keeps the walk sense
//! and is reported as undetermined — the caller never guesses.
//!
//! The loop-sense sign itself is the one every serving ring fixture
//! carries (conic and extreme cube corners, the notched-cube saddle):
//! `(n_F × t) · c < 0`, i.e. a face loop turns clockwise about its
//! outward normal — the kernel-wide loop convention.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, FaceHandle};

use super::boundary_hex::reverse_arc;

/// Ring type produced by `reorient_arcs_for_hex_closure`:
/// `(order, flip, oriented_arcs)`.
pub(crate) type OrientedRing = (Vec<usize>, Vec<bool>, Vec<NurbsCurve>);

/// Outcome of the ring-sense normalization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RingSense {
    /// The walk already circulated in the face-loop sense.
    AlreadyLoopSense,
    /// The walk circulated against it; the cycle was reversed.
    Reversed,
    /// No unanimous bridge vote; walk sense kept.
    Undetermined,
}

/// Per-bridge sense votes: `(loop_sense, mirrored)` counts.
fn bridge_votes(
    brep: &BRep,
    faces: &[FaceHandle],
    vertex_position: Point3,
    oriented_arcs: &[NurbsCurve],
    on_surface_tol: f64,
) -> (usize, usize) {
    let n = oriented_arcs.len();
    let (mut loop_sense, mut mirrored) = (0usize, 0usize);
    for k in 0..n {
        let (_, t1) = oriented_arcs[k].parameter_range();
        let (s0, _) = oriented_arcs[(k + 1) % n].parameter_range();
        let from = oriented_arcs[k].point_at(t1);
        let to = oriented_arcs[(k + 1) % n].point_at(s0);
        let mid = Point3::new(
            0.5 * (from.x + to.x),
            0.5 * (from.y + to.y),
            0.5 * (from.z + to.z),
        );
        let bridge = to - from;
        let Some(vote) = faces.iter().find_map(|&face| {
            let sof = brep.faces.get(face)?.surface.as_ref()?;
            let on = |p: &Point3| {
                sof.geometry.closest_point(p).is_some_and(|(u, v)| {
                    sof.geometry.point_at(u, v).distance_to(p) <= on_surface_tol
                })
            };
            if !(on(&from) && on(&to)) {
                return None;
            }
            let (u, v) = sof.geometry.closest_point(&mid)?;
            let normal = sof.geometry.normal_at(u, v).normalized()?;
            let outward = if sof.orientation {
                normal
            } else {
                normal * -1.0
            };
            let to_vertex = vertex_position - mid;
            let toward_vertex = to_vertex - outward * to_vertex.dot(&outward);
            let side = outward.cross(&bridge).dot(&toward_vertex);
            (side.abs() > f64::EPSILON).then_some(side < 0.0)
        }) else {
            continue;
        };
        if vote {
            loop_sense += 1;
        } else {
            mirrored += 1;
        }
    }
    (loop_sense, mirrored)
}

/// Put a pushed (boundary-ring) cycle into the face-loop sense of
/// the blend face the ring carve will build from it. `faces` are the
/// support faces incident to the vertex; `on_surface_tol` is the
/// bridge builder's on-support tolerance. Trace topic `boundary-hex`
/// records the vote.
pub(crate) fn orient_pushed_ring(
    brep: &BRep,
    faces: &[FaceHandle],
    vertex_position: Point3,
    ring: Option<OrientedRing>,
    on_surface_tol: f64,
) -> Option<OrientedRing> {
    let ring = ring?;
    let (loop_sense, mirrored) =
        bridge_votes(brep, faces, vertex_position, &ring.2, on_surface_tol);
    let (ring, sense) = orient_ring_by_votes(ring, loop_sense, mirrored);
    if popkern_core::trace::trace_enabled("boundary-hex") {
        eprintln!(
            "[boundary-hex] pushed ring sense: {sense:?} (votes loop={loop_sense} mirrored={mirrored})"
        );
    }
    Some(ring)
}

/// Reverse the cycle (order reversed, every arc flipped) exactly when
/// the bridge vote is unanimously mirrored.
fn orient_ring_by_votes(
    ring: OrientedRing,
    loop_sense: usize,
    mirrored: usize,
) -> (OrientedRing, RingSense) {
    match (loop_sense, mirrored) {
        (l, 0) if l > 0 => (ring, RingSense::AlreadyLoopSense),
        (0, m) if m > 0 => {
            let (order, flip, arcs) = ring;
            let order: Vec<usize> = order.into_iter().rev().collect();
            let flip: Vec<bool> = flip.into_iter().rev().map(|f| !f).collect();
            let arcs: Vec<NurbsCurve> = arcs.iter().rev().map(reverse_arc).collect();
            ((order, flip, arcs), RingSense::Reversed)
        }
        _ => (ring, RingSense::Undetermined),
    }
}
