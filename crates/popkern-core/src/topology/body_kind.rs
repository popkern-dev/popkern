//! What kind of body a `BRep` is: empty, an acorn, a wire, a sheet, a
//! solid, or general.
//!
//! One body type holds everything (PK-1101): a body is regions of
//! space, the shells bounding them and the faces, wire edges and
//! isolated vertices those hold. Its kind is read off that structure,
//! cached nowhere, and every operation declares the kinds it accepts
//! ([`BRep::require_kind`]), refusing others with a structured error
//! instead of computing on a body it does not understand:
//!
//! - **Empty**: no faces, edges or vertices.
//! - **Acorn**: isolated vertices only.
//! - **Wire**: edges on no face, no vertex ending more than two of
//!   them.
//! - **Sheet**: faces holding the void on both sides, no edge with
//!   more than two face uses, one fan of faces at each vertex, no wire
//!   edges or acorns.
//! - **Solid**: faces each with solid on one side and void on the
//!   other, every edge with exactly two face uses, one fan at each
//!   vertex, no wire edges or acorns.
//! - **General**: anything else: non-manifold edges or vertices,
//!   internal faces (solid both sides), sheets in or on a solid, lumps
//!   of different kinds, branching wires, a solid opened mid-edit.
//!
//! An edge is **free** when its radial ring holds one face use (a
//! sheet's open boundary; older producers that store it as two uses of
//! one face in two of its loops, the face not closed across the edge,
//! also read free; two uses in one loop are not free, the face lying on
//! both sides: a cylinder's seam, a keyhole bridge, a slit). A shell
//! with a free edge is open, one without closed; "closed" is
//! topological, and a closed solid is not measured for the volume it
//! bounds (validation checks that).

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use super::{BRep, EdgeHandle, RegionKind, ShellHandle};
use crate::error::CoreError;

/// The kind of body a `BRep` is (see the module docs).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BodyKind {
    /// No faces, edges or vertices.
    Empty,
    /// Isolated vertices only.
    Acorn,
    /// Edges on no face, none branching.
    Wire,
    /// Faces with the void on both sides, manifold.
    Sheet,
    /// Faces each between solid and void, closed and manifold.
    Solid,
    /// Anything else: non-manifold, internal faces, mixed lumps.
    General,
}

impl BodyKind {
    /// A short name for messages ("solid", "sheet").
    pub fn label(self) -> &'static str {
        match self {
            BodyKind::Empty => "empty body",
            BodyKind::Acorn => "acorn",
            BodyKind::Wire => "wire",
            BodyKind::Sheet => "sheet",
            BodyKind::Solid => "solid",
            BodyKind::General => "general body",
        }
    }

    /// "solids", "solids or sheets": a list of kinds for messages.
    pub fn list_label(kinds: &[BodyKind]) -> String {
        let names: Vec<String> = kinds.iter().map(|k| format!("{}s", k.label())).collect();
        names.join(" or ")
    }
}

/// How one edge sits in its body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EdgeUse {
    /// Two faces, or one face closed across it (a seam).
    Bounded,
    /// One face holds both half-edges and is open across it.
    Free,
    /// On no face.
    Wire,
}

impl BRep {
    /// The body's kind (see the module docs).
    pub fn body_kind(&self) -> BodyKind {
        if self.faces.is_empty() {
            return if self.edges.is_empty() {
                if self.vertices.is_empty() {
                    BodyKind::Empty
                } else {
                    BodyKind::Acorn
                }
            } else {
                self.wire_kind()
            };
        }
        let mut solid_void = true;
        let mut void_void = true;
        for face in self.faces.values() {
            let solid = |shell| {
                self.shell_region(shell)
                    .and_then(|r| self.region_kind(r))
                    .map(|k| k == RegionKind::Solid)
            };
            match (solid(face.back), solid(face.front)) {
                (Some(false), Some(false)) => solid_void = false,
                (Some(a), Some(b)) if a != b => void_void = false,
                _ => return BodyKind::General,
            }
        }
        let mut rings_of_two = true;
        for edge in self.edges.keys() {
            match self.edge_use(edge) {
                EdgeUse::Wire => return BodyKind::General,
                EdgeUse::Free => rings_of_two = false,
                EdgeUse::Bounded if self.edge_fin_count(edge) > 2 => return BodyKind::General,
                EdgeUse::Bounded => {}
            }
        }
        let has_acorn = self.shells.values().any(|s| s.acorn.is_some());
        let one_fan_each = self.vertex_uses_all().values().all(|uses| uses.len() <= 1);
        if has_acorn || !one_fan_each {
            return BodyKind::General;
        }
        match (solid_void, void_void, rings_of_two) {
            (true, false, true) => BodyKind::Solid,
            (false, true, _) => BodyKind::Sheet,
            _ => BodyKind::General,
        }
    }

    /// A body of edges and no faces: a wire when every edge is a wire
    /// edge, no vertex ends more than two and none is an acorn.
    fn wire_kind(&self) -> BodyKind {
        if self.shells.values().any(|s| s.acorn.is_some()) {
            return BodyKind::General;
        }
        let mut ends: std::collections::HashMap<super::VertexHandle, usize> =
            std::collections::HashMap::new();
        for edge in self.edges.keys() {
            if self.edge_use(edge) != EdgeUse::Wire {
                return BodyKind::General;
            }
            if let Some((a, b)) = self.edge_vertices(edge) {
                *ends.entry(a).or_default() += 1;
                *ends.entry(b).or_default() += 1;
            }
        }
        if ends.values().all(|&n| n <= 2) {
            BodyKind::Wire
        } else {
            BodyKind::General
        }
    }

    /// The gate at an operation's entry: this body's kind when it is one
    /// of `accepted`.
    ///
    /// # Errors
    ///
    /// [`CoreError::UnsupportedBodyKind`] naming the operation, the kind
    /// found and the kinds accepted.
    pub fn require_kind(
        &self,
        operation: &'static str,
        accepted: &'static [BodyKind],
    ) -> Result<BodyKind, CoreError> {
        let found = self.body_kind();
        if accepted.contains(&found) {
            Ok(found)
        } else {
            Err(CoreError::UnsupportedBodyKind {
                operation,
                found,
                accepted,
            })
        }
    }

    /// The shells with a free edge on one of their faces: a sheet's
    /// shells, open where a solid's would close.
    pub(crate) fn open_shells(&self) -> BTreeSet<ShellHandle> {
        self.edges
            .keys()
            .filter(|&e| self.edge_use(e) == EdgeUse::Free)
            .filter_map(|e| self.edge_shell(e))
            .collect()
    }

    /// Whether `edge` is a free (open boundary) edge: its fin is alone
    /// on it, a ring of one (null ring links). An edge whose fin has a
    /// mate is not free even while the mate is in no loop yet (a fin
    /// being wired, a wire edge's half-edge pair); a seam (two uses of
    /// one periodic face) is not free.
    pub fn is_free_edge(&self, edge: EdgeHandle) -> bool {
        use slotmap::Key;
        self.edges
            .get(edge)
            .and_then(|e| self.half_edges.get(e.half_edge))
            .is_some_and(|fin| fin.ring.next.is_null() && fin.ring.prev.is_null())
    }

    /// Whether the body has any free edge (an open sheet boundary).
    pub fn has_free_edges(&self) -> bool {
        self.edges.keys().any(|e| self.is_free_edge(e))
    }

    fn edge_use(&self, edge: EdgeHandle) -> EdgeUse {
        if self.is_wire_edge(edge) {
            return EdgeUse::Wire;
        }
        let fins = self.radial_fins(edge);
        match fins.as_slice() {
            // No face uses and not a wire edge: a broken edge (a missing
            // half-edge reads as free so a broken body never reads as
            // closed).
            [] | [_] => EdgeUse::Free,
            // Two uses of one face in two of its loops, the face not
            // closing across the edge: the older mirrored sheet
            // boundary (a face's own loop and a mirror of it). Two uses
            // in ONE loop have the face on both sides: a bridge the
            // loop crosses twice (a keyhole joining an annulus's rims,
            // a seam) or a spur the loop runs out and back along (a
            // slit a sheet leaves dangling in a solid's face, mvfs +
            // mev's edge), neither of which opens the shell.
            [a, _] if self.is_mirrored_boundary_edge(edge) => {
                let face = self
                    .half_edges
                    .get(*a)
                    .and_then(|h| self.loops.get(h.loop_))
                    .map(|l| l.face);
                if face.is_some_and(|f| !self.face_closed_across_seam(f)) {
                    EdgeUse::Free
                } else {
                    EdgeUse::Bounded
                }
            }
            _ => EdgeUse::Bounded,
        }
    }

    fn edge_shell(&self, edge: EdgeHandle) -> Option<ShellHandle> {
        let e = self.edges.get(edge)?;
        let h = self.half_edges.get(e.half_edge)?;
        let l = self.loops.get(h.loop_)?;
        self.faces.get(l.face).map(|f| f.back)
    }

    /// Whether `face`'s surface closes on itself, so an edge holding the
    /// face on both sides is its seam: the analytic periodic surfaces,
    /// and a NURBS surface whose opposite parameter edges coincide (a
    /// swept circular profile's side face).
    fn face_closed_across_seam(&self, face: super::FaceHandle) -> bool {
        use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
        let Some(surf) = self.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            return false;
        };
        match &surf.geometry {
            SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Cone { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Torus { .. } => true,
            SurfaceGeometry::Nurbs(_) => {
                let s = &surf.geometry;
                let ((u0, u1), (v0, v1)) = s.parameter_range();
                let gate = self.tolerance.linear().max(1.0e-9) * 10.0;
                let closed_u = (0..=4).all(|i| {
                    let v = v0 + (v1 - v0) * i as f64 / 4.0;
                    s.point_at(u0, v).distance_to(&s.point_at(u1, v)) < gate
                });
                let closed_v = (0..=4).all(|i| {
                    let u = u0 + (u1 - u0) * i as f64 / 4.0;
                    s.point_at(u, v0).distance_to(&s.point_at(u, v1)) < gate
                });
                closed_u || closed_v
            }
            _ => false,
        }
    }
}

#[cfg(test)]
#[path = "body_kind_tests.rs"]
mod tests;
