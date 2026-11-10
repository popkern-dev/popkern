//! A face the draft consumes, served as the ridge where the drafted
//! faces beside it meet (PK-1159, the tangency foundations' draft
//! consumer, design §8.4, F1).
//!
//! Drafting both long sides of a rib towards each other until they
//! meet below its top turns the top face's loop inside out: its
//! recomputed vertices lie past each other. The body the draft means
//! has no top face; the sides meet in a ridge, and each end face's top
//! edge shrinks to the apex where the ridge meets it. That topology is
//! made on the input before the draft runs again: the strip face's two
//! short edges are contracted (`kve`) and the face is merged into one
//! drafted side (`kef`), leaving its other long edge as the ridge
//! between the two sides. The draft then places the apexes as the
//! meeting of the two drafted surfaces with the end face, and the
//! ridge as the drafted surfaces' meeting, as it places every vertex
//! and edge.
//!
//! The class served: a planar strip of four edges, two opposite ones
//! shared with two distinct drafted faces. Any other consumed face
//! still declines `FaceConsumed`.

use std::collections::HashMap;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// `brep` with the consumed strip `face` collapsed into a ridge between
/// the drafted faces across its long edges, or `None` when the face is
/// not such a strip or an operator declines.
pub(crate) fn collapse_consumed_strip(
    brep: &BRep,
    face: FaceHandle,
    drafted: &HashMap<FaceHandle, SurfaceGeometry>,
) -> Option<BRep> {
    let f = brep.faces.get(face)?;
    if !f.inner_loops.is_empty() {
        return None;
    }
    let start = brep.loops.get(f.outer_loop)?.half_edge;
    let mut ring: Vec<HalfEdgeHandle> = Vec::new();
    let mut he = start;
    for _ in 0..=4 {
        ring.push(he);
        he = brep.half_edges.get(he)?.next;
        if he == start {
            break;
        }
    }
    if ring.len() != 4 || he != start {
        return None;
    }
    let across = |h: HalfEdgeHandle| -> Option<FaceHandle> {
        let mate = brep.half_edges.get(h)?.mate()?;
        Some(brep.loops.get(brep.half_edges.get(mate)?.loop_)?.face)
    };
    let sides: Vec<Option<FaceHandle>> = ring.iter().map(|&h| across(h)).collect();
    let drafted_side = |k: usize| sides[k].filter(|s| drafted.contains_key(s) && *s != face);
    let (long, short) = [(0usize, 1usize), (1, 0)].into_iter().find(
        |&(l, _)| matches!((drafted_side(l), drafted_side(l + 2)), (Some(a), Some(b)) if a != b),
    )?;
    let (long1, short1, short2) = (ring[long], ring[short], ring[short + 2]);
    let mut work = brep.clone();
    work.kve(short1).ok()?;
    // The strip merges into the side across its first long edge: kef
    // keeps the half-edge's face and kills its twin's.
    let side_half = work.half_edges.get(long1)?.mate()?;
    work.kef(side_half).ok()?;
    work.kve(short2).ok()?;
    Some(work)
}
