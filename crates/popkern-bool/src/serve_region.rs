//! The serve tail's touched region (PK-815,
//! docs/architecture/local-boolean-scope.md "the touched region").
//!
//! A locally scoped boolean carries most of `a` through unchanged. Its
//! serve tail normalizes and re-derives only the TOUCHED faces: those
//! with an edge not built from an unsplit original of `a`. Structure
//! carried unchanged keeps its form and derived state (the locality
//! decision). Validity is not regional: a pass that edits drops the
//! kept verdict and the next gate re-validates the whole body. A
//! result that carried nothing (unscoped, graph-merge) is touched
//! everywhere, and every pass runs on the whole body as before.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, LoopHandle, VertexHandle};
use slotmap::SecondaryMap;

/// Result edges built from unsplit original edges. Dense on the
/// result's handles: the serve tail tests every edge.
pub(crate) type CarriedEdges = SecondaryMap<EdgeHandle, ()>;

/// What a cellular result carries over from its operands.
#[derive(Default)]
pub(crate) struct Carried {
    /// Edges built from unsplit original edges of either operand: their
    /// curve and both faces' surfaces are an operand's own, so the
    /// curve-surface check and the tolerance measurement skip them.
    pub(crate) edges: CarriedEdges,
    /// The edges built from `a`'s unsplit edges.
    pub(crate) a_edges: CarriedEdges,
    /// The faces built from `a`'s (kept fragments, carried far faces).
    pub(crate) a_faces: SecondaryMap<FaceHandle, ()>,
    /// The edges built from `b`'s unsplit edges, recorded only when its
    /// far faces are carried (a scoped union, PK-816).
    pub(crate) b_edges: CarriedEdges,
    /// The faces built from `b`'s kept fragments and carried far faces,
    /// under the same condition, when they stand as in `b` (unflipped,
    /// or all reversed together to agree with `a`, PK-916).
    pub(crate) b_faces: SecondaryMap<FaceHandle, ()>,
}

impl Carried {
    /// One operand's carried edges and faces.
    pub(crate) fn of(
        &mut self,
        operand: crate::merge::graph::Operand,
    ) -> (&mut CarriedEdges, &mut SecondaryMap<FaceHandle, ()>) {
        match operand {
            crate::merge::graph::Operand::A => (&mut self.a_edges, &mut self.a_faces),
            crate::merge::graph::Operand::B => (&mut self.b_edges, &mut self.b_faces),
        }
    }
}

/// Where the serve tail may edit.
pub(crate) enum ServeRegion {
    /// Every face (nothing was carried).
    Whole,
    /// The touched faces and loops.
    Touched(Touched),
}

/// The faces and loops holding an uncarried edge.
pub(crate) struct Touched {
    /// Face membership.
    faces: SecondaryMap<FaceHandle, ()>,
    /// The same faces, in handle order.
    list: Vec<FaceHandle>,
    /// Loop membership.
    loops: SecondaryMap<LoopHandle, ()>,
}

impl ServeRegion {
    /// The region of `brep` given what it `carried` from `a`.
    pub(crate) fn of(brep: &BRep, carried: &Carried) -> Self {
        if carried.a_edges.is_empty() && carried.b_edges.is_empty() {
            return Self::Whole;
        }
        Self::marked(brep, carried, SecondaryMap::new(), SecondaryMap::new())
    }

    /// Re-derive after an edit: faces and loops created or newly
    /// holding an uncarried edge join; those already in keep their
    /// place (an edit can leave what it changed with carried edges
    /// only), unless gone.
    pub(crate) fn refresh(&mut self, brep: &BRep, carried: &Carried) {
        let Self::Touched(t) = self else {
            return;
        };
        let faces = std::mem::take(&mut t.faces)
            .into_iter()
            .filter(|&(f, ())| brep.faces.contains_key(f))
            .collect();
        let loops = std::mem::take(&mut t.loops)
            .into_iter()
            .filter(|&(l, ())| brep.loops.contains_key(l))
            .collect();
        *self = Self::marked(brep, carried, faces, loops);
    }

    /// Add both loops and faces of every edge that is not INTACT: an
    /// edge is intact when it is built from an unsplit edge of `a` and
    /// both its faces are built from `a`'s, i.e. it and its two sides
    /// stand as in `a` — or the same of `b` when a scoped union carries
    /// `b`'s far faces (PK-816). A face or loop of intact edges only
    /// keeps an operand's structure (a carried hole, an unsplit boundary
    /// between unsplit faces); a face along the operands' meeting, a
    /// flush neighbour that replaced one of `a`'s, and every face a pass
    /// creates have an edge that is not.
    fn marked(
        brep: &BRep,
        carried: &Carried,
        mut faces: SecondaryMap<FaceHandle, ()>,
        mut loops: SecondaryMap<LoopHandle, ()>,
    ) -> Self {
        for (eh, edge) in &brep.edges {
            let Some(he) = brep.half_edges.get(edge.half_edge) else {
                continue;
            };
            let sides = [he, brep.half_edges.get(he.twin).unwrap_or(he)];
            let face_of =
                |h: &popkern_core::topology::HalfEdge| brep.loops.get(h.loop_).map(|l| l.face);
            let stands = |edges: &CarriedEdges, faces: &SecondaryMap<FaceHandle, ()>| {
                edges.contains_key(eh)
                    && sides
                        .iter()
                        .all(|h| face_of(h).is_some_and(|f| faces.contains_key(f)))
            };
            let intact = stands(&carried.a_edges, &carried.a_faces)
                || stands(&carried.b_edges, &carried.b_faces);
            if intact {
                continue;
            }
            for h in sides {
                if let Some(f) = face_of(h) {
                    loops.insert(h.loop_, ());
                    faces.insert(f, ());
                }
            }
        }
        let list = faces.keys().collect();
        Self::Touched(Touched { faces, list, loops })
    }

    /// Whether `face` may be edited.
    pub(crate) fn admits_face(&self, face: FaceHandle) -> bool {
        match self {
            Self::Whole => true,
            Self::Touched(t) => t.faces.contains_key(face),
        }
    }

    /// Whether `l` can hold a candidate a loop-local pass would edit:
    /// a loop of carried edges only keeps `a`'s structure there.
    pub(crate) fn admits_loop(&self, l: LoopHandle) -> bool {
        match self {
            Self::Whole => true,
            Self::Touched(t) => t.loops.contains_key(l),
        }
    }

    /// Whether an edit of `edge` (its two faces) stays in the region.
    pub(crate) fn admits_edge(&self, brep: &BRep, edge: EdgeHandle) -> bool {
        if matches!(self, Self::Whole) {
            return true;
        }
        let Some(he) = brep
            .edges
            .get(edge)
            .and_then(|e| brep.half_edges.get(e.half_edge))
        else {
            return false;
        };
        [Some(he), brep.half_edges.get(he.twin)]
            .into_iter()
            .all(|h| {
                h.and_then(|h| brep.loops.get(h.loop_))
                    .is_some_and(|l| self.admits_face(l.face))
            })
    }

    /// Whether an edit at `vertex` (every face a half-edge leaves it
    /// in) stays in the region. A scan of every half-edge, exact on any
    /// wiring (a fan walk breaks at the self-twins the weld passes
    /// repair): for candidate sites, not per-vertex screens.
    pub(crate) fn admits_vertex(&self, brep: &BRep, vertex: VertexHandle) -> bool {
        matches!(self, Self::Whole)
            || brep
                .half_edges
                .values()
                .filter(|h| h.origin == vertex)
                .all(|h| {
                    brep.loops
                        .get(h.loop_)
                        .is_some_and(|l| self.admits_face(l.face))
                })
    }

    /// Re-derive the region's orientation flags; others keep theirs.
    pub(crate) fn reconcile_orientation_flags(&self, brep: &mut BRep) -> usize {
        match self {
            Self::Whole => brep.reconcile_orientation_flags(),
            Self::Touched(t) => brep.reconcile_orientation_flags_among(&t.list),
        }
    }
}
