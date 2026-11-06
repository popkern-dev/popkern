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
//! The cap builders on the ring construct each sector's patch from the
//! ring's circulation, so the ring must turn the way its builder needs
//! about the MATERIAL: an exact local reading at every bridge rather
//! than an averaged direction (a saddle's support normals sum to the
//! concave side and mislead any average). Bridge `k` runs across
//! support face `F` from arc `k`'s end to arc `k+1`'s start, and the
//! carve hands the vertex side of that cut to the cap. With `m_F` the
//! direction into the material at `F` (its flag-folded normal turned
//! by the body's material sense), `t` the bridge direction and `c` the
//! in-surface direction from the bridge toward the vertex,
//! `(m_F × t) · c < 0` for every bridge of a ring turning
//! counter-clockwise about the into-material direction and `> 0` for
//! every bridge of one turning clockwise. The bridges vote; a unanimous
//! vote decides. A split or empty vote (no support face found,
//! degenerate bridge, a material sense that cannot be measured) keeps
//! the walk sense and is reported as undetermined: the caller never
//! guesses.
//!
//! The convex cap builders (N-quad Gregory, PnS) need the ring
//! counter-clockwise about the into-material direction, the saddle
//! N-quad cap clockwise (PK-523: both measured on bodies whose
//! flag-folded normals pointed into the material, where the vote read
//! the folded normal; the saddle's cap near-folds in the other sense,
//! in-patch normal cos -0.993 against the -0.98 reversal gate). Read
//! against the material, the same corner of a body facing either way
//! gets the same ring.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, FaceHandle};

use super::boundary_hex::reverse_arc;

/// Ring type produced by `reorient_arcs_for_hex_closure`:
/// `(order, flip, oriented_arcs)`.
pub(crate) type OrientedRing = (Vec<usize>, Vec<bool>, Vec<NurbsCurve>);

/// The sense a ring's cap builder needs, about the direction into the
/// material (module docs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingSenseRule {
    /// Counter-clockwise about the into-material direction: the convex
    /// ring rung's cap builders.
    Convex,
    /// Clockwise about it: the saddle N-quad cap.
    Saddle,
}

/// Outcome of the ring-sense normalization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RingSense {
    /// The walk already circulated in the rule's sense.
    AlreadyInSense,
    /// The walk circulated against it; the cycle was reversed.
    Reversed,
    /// No unanimous bridge vote; walk sense kept.
    Undetermined,
}

/// Per-bridge sense votes: `(counter_clockwise, clockwise)` counts
/// about the into-material direction, `material_out` the body's
/// out-of-material multiplier for flag-folded normals.
fn bridge_votes(
    brep: &BRep,
    faces: &[FaceHandle],
    vertex_position: Point3,
    oriented_arcs: &[NurbsCurve],
    on_surface_tol: f64,
    material_out: f64,
) -> (usize, usize) {
    let n = oriented_arcs.len();
    let (mut counter_clockwise, mut clockwise) = (0usize, 0usize);
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
            let folded = if sof.orientation {
                normal
            } else {
                normal * -1.0
            };
            let into = folded * -material_out;
            let to_vertex = vertex_position - mid;
            let toward_vertex = to_vertex - into * to_vertex.dot(&into);
            let side = into.cross(&bridge).dot(&toward_vertex);
            (side.abs() > f64::EPSILON).then_some(side < 0.0)
        }) else {
            continue;
        };
        if vote {
            counter_clockwise += 1;
        } else {
            clockwise += 1;
        }
    }
    (counter_clockwise, clockwise)
}

/// Put a boundary-ring cycle into the sense `rule` asks of it about
/// the material. `faces` are the support faces incident to the vertex;
/// `on_surface_tol` is the bridge builder's on-support tolerance. Trace
/// topic `boundary-hex` records the vote.
pub(crate) fn orient_pushed_ring(
    brep: &BRep,
    faces: &[FaceHandle],
    vertex_position: Point3,
    ring: Option<OrientedRing>,
    on_surface_tol: f64,
    rule: RingSenseRule,
) -> Option<OrientedRing> {
    let ring = ring?;
    let material_out = faces
        .first()
        .and_then(|&f| brep.component_material_sense(f))
        .map(|s| s.outward_multiplier());
    let (ccw, cw) = match material_out {
        Some(m) => bridge_votes(brep, faces, vertex_position, &ring.2, on_surface_tol, m),
        None => (0, 0),
    };
    let (in_sense, against) = match rule {
        RingSenseRule::Convex => (ccw, cw),
        RingSenseRule::Saddle => (cw, ccw),
    };
    let (ring, sense) = orient_ring_by_votes(ring, in_sense, against);
    if popkern_core::trace::trace_enabled("boundary-hex") {
        eprintln!(
            "[boundary-hex] ring sense {rule:?}: {sense:?} (votes counter-clockwise={ccw} clockwise={cw}, material {material_out:?})"
        );
    }
    Some(ring)
}

/// Reverse the cycle (order reversed, every arc flipped) exactly when
/// the bridge vote is unanimously against the wanted sense.
fn orient_ring_by_votes(
    ring: OrientedRing,
    in_sense: usize,
    against: usize,
) -> (OrientedRing, RingSense) {
    match (in_sense, against) {
        (l, 0) if l > 0 => (ring, RingSense::AlreadyInSense),
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
