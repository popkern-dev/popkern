//! A closed rim that lies about its direction must not convict its
//! honest neighbours (PK-741, from PK-730's finned housing: a
//! seam-fin union's top rim curve belonged to the cap, the wall's half-
//! edge was made canonical, the wall's loop wound 4π, the periodic
//! closure flipped the honest bottom rim, and the canonical-forward
//! audit then convicted the honest cap).
//!
//! The plain cylinder carries the whole structure: a periodic wall
//! loop with two full-period rims, each shared with a planar cap.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};
use popkern_core::validation::CanonicalForwardState;

fn cylinder() -> BRep {
    BRep::cylinder(2.0, 3.0, Tolerance::default()).expect("cylinder")
}

fn wall(b: &BRep) -> FaceHandle {
    b.faces
        .iter()
        .find(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| s.geometry.kind_name() == "Cylinder")
        })
        .map(|(h, _)| h)
        .expect("wall")
}

/// The closed rims, with each one's wall-side half-edge.
fn rims(b: &BRep) -> Vec<(EdgeHandle, HalfEdgeHandle)> {
    let wall = wall(b);
    b.edges
        .iter()
        .filter_map(|(eh, e)| {
            let he = &b.half_edges[e.half_edge];
            if he.origin != b.half_edges[he.twin].origin {
                return None;
            }
            let on_wall = |h: HalfEdgeHandle| b.loops[b.half_edges[h].loop_].face == wall;
            let wall_he = [e.half_edge, he.twin].into_iter().find(|&h| on_wall(h))?;
            Some((eh, wall_he))
        })
        .collect()
}

fn top_rim(b: &BRep) -> (EdgeHandle, HalfEdgeHandle) {
    use popkern_core::geometry::CurveEval;
    rims(b)
        .into_iter()
        .find(|&(eh, _)| {
            b.edges[eh]
                .curve
                .as_ref()
                .is_some_and(|c| c.geometry.point_at(c.t_start).z > 1.5)
        })
        .expect("top rim")
}

/// The wall-side rim pcurves' running direction in u, per rim.
fn wall_rim_directions(b: &BRep) -> Vec<(EdgeHandle, bool)> {
    rims(b)
        .into_iter()
        .map(|(eh, wall_he)| {
            let pc = b.half_edges[wall_he].pcurve.as_ref().expect("pcurve");
            let (t0, t1) = pc.parameter_range();
            (eh, pc.point_at(t1).x > pc.point_at(t0).x)
        })
        .collect()
}

fn assert_honest(b: &BRep, clean: &BRep) {
    for v in b.audit_canonical_forward() {
        assert_eq!(v.state, CanonicalForwardState::Ok, "edge {:?}", v.edge);
    }
    assert!(
        b.dishonest_winding_faces().is_empty(),
        "winding-dishonest faces"
    );
    assert_eq!(wall_rim_directions(b), wall_rim_directions(clean));
    for (fh, f) in &b.faces {
        assert_eq!(
            f.surface.as_ref().map(|s| s.orientation),
            clean.faces[fh].surface.as_ref().map(|s| s.orientation),
            "face {fh:?} flag"
        );
    }
}

/// The PK-730 lie: the top rim's wall half-edge made canonical with
/// its pcurve claiming the cap's direction, so both sides claim
/// forward and the wall loop winds 4π. Recompute and repair restore
/// the honest rim. The closure must flip the LYING member, the one
/// whose claim agrees with its twin's, not the first full-period one.
#[test]
fn a_lying_rim_is_flipped_and_repaired_not_its_honest_sibling() {
    let clean = cylinder();
    let mut b = clean.clone();
    let (top, wall_he) = top_rim(&b);
    b.edges[top].half_edge = wall_he;
    let lie = b.half_edges[wall_he]
        .pcurve
        .as_ref()
        .expect("pcurve")
        .reversed();
    b.half_edges[wall_he].pcurve = Some(lie);
    let wall = wall(&b);
    b.recompute_face_pcurves(wall).expect("recompute");
    // The bottom rim is untouched by the closure.
    let bottom_dirs: Vec<_> = wall_rim_directions(&b)
        .into_iter()
        .filter(|&(eh, _)| eh != top)
        .collect();
    let clean_bottom: Vec<_> = wall_rim_directions(&clean)
        .into_iter()
        .filter(|&(eh, _)| eh != top)
        .collect();
    assert_eq!(
        bottom_dirs, clean_bottom,
        "the honest bottom rim was flipped"
    );
    b.repair_canonical_forward();
    assert_honest(&b, &clean);
}

/// A wall whose rims are BOTH reversed and whose seams were then
/// re-authored to match is self-consistent: the chain test passes it.
/// Its loop winds against its own flag, though, while the caps agree
/// with theirs, so the audit convicts the wall's pcurves on both rims
/// and never a cap's.
#[test]
fn a_self_consistent_reversed_wall_is_the_convicted_side() {
    let mut b = cylinder();
    for (eh, wall_he) in rims(&b) {
        let reversed = b.half_edges[wall_he]
            .pcurve
            .as_ref()
            .expect("pcurve")
            .reversed();
        b.half_edges[wall_he].pcurve = Some(reversed);
        b.edges[eh].half_edge = wall_he;
    }
    let wall = wall(&b);
    b.recompute_face_pcurves(wall).expect("recompute");
    let walls: Vec<HalfEdgeHandle> = rims(&b).into_iter().map(|(_, h)| h).collect();
    for v in b.audit_canonical_forward() {
        let liar = match v.state {
            CanonicalForwardState::OkByChain { lying_half_edge }
            | CanonicalForwardState::ViolatedByChain { lying_half_edge } => lying_half_edge,
            other => panic!("edge {:?}: {other:?}", v.edge),
        };
        assert!(walls.contains(&liar), "edge {:?} convicted the cap", v.edge);
    }
}
