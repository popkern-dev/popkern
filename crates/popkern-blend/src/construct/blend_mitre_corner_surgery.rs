//! Euler surgery for the mitre-corner settlement.
//!
//! Two phases, both through Euler operators only:
//!
//! 1. [`reduce_to_residue_polygon`]: the flaps around the corner
//!    (every face incident to it) merge into one residue face (`kef`
//!    across flap|flap edges), and the spokes left dangling die (`kev`).
//!    The residue face must then be bounded only by strip edges, each
//!    strip in one contiguous run. The vertices where the bordering
//!    strip changes are the rail crossings. Each run's interior
//!    scaffold vertices contract into its first crossing (`kve`),
//!    leaving an N-gon with one strip across each side. The residue
//!    shape does not matter. Orthogonal corners leave quad flaps and
//!    oblique ones leave chains of stub vertices; both reduce the same
//!    way.
//! 2. [`star_from_mitre_point`]: a vertex at the mitre point joins
//!    every crossing (`mev`, then `mef` per remaining crossing). Each
//!    triangle then merges into the strip across its polygon side
//!    (`kef`), which ends every strip in a V at the mitre point.
//!
//! The caller owns rollback: a failed phase leaves a half-built body
//! that the caller restores from its snapshot.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use crate::types::MitreCornerDeclineReason as Reason;

/// The reduced corner: an N-gon of rail crossings.
pub(super) struct ResiduePolygon {
    /// The residue face.
    pub(super) face: FaceHandle,
    /// Rail crossings in the residue loop's order.
    pub(super) crossings: Vec<VertexHandle>,
    /// `sides[i]` is the strip across the side `crossings[i] →
    /// crossings[i + 1]`.
    pub(super) sides: Vec<FaceHandle>,
}

fn face_of(brep: &BRep, he: HalfEdgeHandle) -> Option<FaceHandle> {
    Some(brep.loops.get(brep.half_edges.get(he)?.loop_)?.face)
}

fn loop_ring(brep: &BRep, loop_: LoopHandle) -> Option<Vec<HalfEdgeHandle>> {
    let start = brep.loops.get(loop_)?.half_edge;
    let mut ring = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        ring.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(ring);
        }
    }
    None
}

fn outer_ring(brep: &BRep, face: FaceHandle) -> Option<Vec<HalfEdgeHandle>> {
    let face_data = brep.faces.get(face)?;
    if !face_data.inner_loops.is_empty() {
        return None;
    }
    loop_ring(brep, face_data.outer_loop)
}

/// Phase 1 (see the module docs).
pub(super) fn reduce_to_residue_polygon(
    brep: &mut BRep,
    corner: VertexHandle,
    strips: &[FaceHandle],
) -> Result<ResiduePolygon, Reason> {
    let mut flaps: Vec<FaceHandle> = Vec::new();
    for (h, he) in brep.half_edges.iter() {
        if he.origin != corner {
            continue;
        }
        let face = face_of(brep, h).ok_or(Reason::UnrecognizedResidue)?;
        if strips.contains(&face) {
            return Err(Reason::UnrecognizedResidue);
        }
        if !flaps.contains(&face) {
            flaps.push(face);
        }
    }
    // Merge flaps across their shared edges until one face remains.
    loop {
        let shared = brep.half_edges.iter().find_map(|(h, he)| {
            let own = face_of(brep, h)?;
            let across = face_of(brep, he.twin)?;
            (own != across && flaps.contains(&own) && flaps.contains(&across)).then_some(h)
        });
        let Some(h) = shared else {
            break;
        };
        let killed = face_of(brep, brep.half_edges[h].twin).ok_or(Reason::EulerRejected)?;
        brep.kef(h).map_err(|_| Reason::EulerRejected)?;
        flaps.retain(|f| *f != killed);
    }
    let [residue] = flaps[..] else {
        return Err(Reason::UnrecognizedResidue);
    };
    // Retire dangling spokes tip-first: an edge with the residue on
    // both sides whose far vertex has no other edge.
    loop {
        let ring = outer_ring(brep, residue).ok_or(Reason::UnrecognizedResidue)?;
        let spur = ring.iter().copied().find(|&h| {
            let he = &brep.half_edges[h];
            he.next == he.twin
        });
        let Some(h) = spur else {
            break;
        };
        brep.kev(h).map_err(|_| Reason::EulerRejected)?;
    }

    // The residue must now be bounded by strip runs only.
    let ring = outer_ring(brep, residue).ok_or(Reason::UnrecognizedResidue)?;
    let across: Vec<FaceHandle> = ring
        .iter()
        .map(|&h| face_of(brep, brep.half_edges[h].twin).ok_or(Reason::UnrecognizedResidue))
        .collect::<Result<_, _>>()?;
    if across.iter().any(|f| !strips.contains(f)) {
        return Err(Reason::UnrecognizedResidue);
    }
    // Run starts: positions where the strip across changes.
    let n = ring.len();
    let starts: Vec<usize> = (0..n)
        .filter(|&i| across[i] != across[(i + n - 1) % n])
        .collect();
    let mut sides: Vec<FaceHandle> = starts.iter().map(|&i| across[i]).collect();
    let mut sorted = sides.clone();
    sorted.sort();
    sorted.dedup();
    if starts.len() != strips.len() || sorted.len() != sides.len() {
        return Err(Reason::UnrecognizedResidue);
    }
    let crossings: Vec<VertexHandle> = starts
        .iter()
        .map(|&i| brep.half_edges[ring[i]].origin)
        .collect();

    // Contract every run to one side: merge each interior vertex into
    // the run's first crossing (kve keeps the origin's position).
    for &x in &crossings {
        loop {
            let ring = outer_ring(brep, residue).ok_or(Reason::EulerRejected)?;
            let Some(&h) = ring.iter().find(|&&h| brep.half_edges[h].origin == x) else {
                return Err(Reason::EulerRejected);
            };
            let dest = brep.half_edges[brep.half_edges[h].twin].origin;
            if crossings.contains(&dest) {
                break;
            }
            brep.kve(h).map_err(|_| Reason::EulerRejected)?;
        }
    }
    let ring = outer_ring(brep, residue).ok_or(Reason::EulerRejected)?;
    if ring.len() != crossings.len() {
        return Err(Reason::UnrecognizedResidue);
    }
    // Re-read the order and sides from the reduced loop.
    let start = ring
        .iter()
        .position(|&h| brep.half_edges[h].origin == crossings[0])
        .ok_or(Reason::EulerRejected)?;
    let mut ordered = Vec::with_capacity(ring.len());
    sides.clear();
    for k in 0..ring.len() {
        let h = ring[(start + k) % ring.len()];
        ordered.push(brep.half_edges[h].origin);
        sides.push(face_of(brep, brep.half_edges[h].twin).ok_or(Reason::EulerRejected)?);
    }
    Ok(ResiduePolygon {
        face: residue,
        crossings: ordered,
        sides,
    })
}

/// Phase 2 (see the module docs). Returns the mitre vertex.
pub(super) fn star_from_mitre_point(
    brep: &mut BRep,
    polygon: &ResiduePolygon,
    mitre_point: Point3,
) -> Result<VertexHandle, Reason> {
    let side_half_edge = |brep: &BRep, from: VertexHandle, strip: FaceHandle| {
        brep.half_edges.iter().find_map(|(h, he)| {
            (he.origin == from && face_of(brep, he.twin) == Some(strip)).then_some(h)
        })
    };
    let first = side_half_edge(brep, polygon.crossings[0], polygon.sides[0])
        .ok_or(Reason::EulerRejected)?;
    let spur = brep
        .mev(first, mitre_point)
        .map_err(|_| Reason::EulerRejected)?;
    let mitre = spur.vertex;
    for (k, &x) in polygon.crossings.iter().enumerate().skip(1) {
        let to = side_half_edge(brep, x, polygon.sides[k]).ok_or(Reason::EulerRejected)?;
        let loop_ = brep.half_edges[to].loop_;
        let from = brep
            .half_edges
            .iter()
            .find_map(|(h, he)| (he.origin == mitre && he.loop_ == loop_).then_some(h))
            .ok_or(Reason::EulerRejected)?;
        brep.mef(from, to).map_err(|_| Reason::EulerRejected)?;
    }
    // Every triangle joins the strip across its polygon side.
    for (k, &x) in polygon.crossings.iter().enumerate() {
        let side = side_half_edge(brep, x, polygon.sides[k]).ok_or(Reason::EulerRejected)?;
        let strip_side = brep.half_edges[side].twin;
        brep.kef(strip_side).map_err(|_| Reason::EulerRejected)?;
    }
    if brep.faces.get(polygon.face).is_some() {
        return Err(Reason::EulerRejected);
    }
    Ok(mitre)
}
