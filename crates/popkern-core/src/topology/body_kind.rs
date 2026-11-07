//! What kind of body a `BRep` is: solid, sheet, wire, a mix, or empty.
//!
//! Commercial kernels type their bodies and every operation declares
//! the kinds it accepts (plan D4): a boolean of solids refuses a sheet
//! with a structured error instead of classifying volume it does not
//! bound. The kind is decided from topology alone, per lump (shell):
//!
//! - an edge is **free** when its radial ring holds one face use (a
//!   sheet's open boundary; legacy producers that store it as two uses
//!   of one face in two of its loops, the face not closed across the
//!   edge, also read free; two uses in one loop are not free, the face
//!   lying on both sides: a cylinder's seam, a keyhole bridge, a slit);
//! - an edge on no face is a **wire** edge;
//! - a shell with a free edge is a sheet, one without is closed.
//!
//! "Closed" is topological: a closed shell is reported `Solid` without
//! measuring the volume it bounds.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use super::{BRep, EdgeHandle, ShellHandle};

/// The kind of body a `BRep` is.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BodyKind {
    /// No faces and no edges.
    Empty,
    /// Every shell closed.
    Solid,
    /// Faces only, at least one shell open (free edges).
    Sheet,
    /// Edges on no face, and no faces.
    Wire,
    /// More than one of the kinds above among its lumps.
    Mixed,
}

impl BodyKind {
    /// A short name for messages ("solid", "sheet").
    pub fn label(self) -> &'static str {
        match self {
            BodyKind::Empty => "empty body",
            BodyKind::Solid => "solid",
            BodyKind::Sheet => "sheet",
            BodyKind::Wire => "wire",
            BodyKind::Mixed => "mixed body",
        }
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
                BodyKind::Empty
            } else {
                BodyKind::Wire
            };
        }
        let wire = self.edges.keys().any(|e| self.edge_use(e) == EdgeUse::Wire);
        let open = self.open_shells();
        let shells: BTreeSet<ShellHandle> = self.faces.values().map(|f| f.back).collect();
        let sheet = !open.is_empty();
        let solid = shells.iter().any(|s| !open.contains(s));
        match (solid, sheet, wire) {
            (true, false, false) => BodyKind::Solid,
            (false, true, false) => BodyKind::Sheet,
            _ => BodyKind::Mixed,
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
            // closing across the edge: the legacy mirrored sheet
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
