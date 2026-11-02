//! What kind of body a `BRep` is: solid, sheet, wire, a mix, or empty.
//!
//! Commercial kernels type their bodies and every operation declares
//! the kinds it accepts (plan D4): a boolean of solids refuses a sheet
//! with a structured error instead of classifying volume it does not
//! bound. The kind is decided from topology alone, per lump (shell):
//!
//! - an edge is **free** when one face holds both its half-edges and
//!   that face's surface is not closed across the edge (a sheet's open
//!   boundary is stored that way; a cylinder's seam is not free);
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
        let mut open: BTreeSet<ShellHandle> = BTreeSet::new();
        let mut wire = false;
        for (eh, _) in &self.edges {
            match self.edge_use(eh) {
                EdgeUse::Bounded => {}
                EdgeUse::Wire => wire = true,
                EdgeUse::Free => {
                    if let Some(shell) = self.edge_shell(eh) {
                        open.insert(shell);
                    }
                }
            }
        }
        let shells: BTreeSet<ShellHandle> = self.faces.values().map(|f| f.shell).collect();
        let sheet = !open.is_empty();
        let solid = shells.iter().any(|s| !open.contains(s));
        match (solid, sheet, wire) {
            (true, false, false) => BodyKind::Solid,
            (false, true, false) => BodyKind::Sheet,
            _ => BodyKind::Mixed,
        }
    }

    /// Whether `edge` is a free (open boundary) edge: one face holds
    /// both its half-edges and that face's surface is not closed across
    /// it. An edge with a missing half-edge, twin or loop counts as free
    /// when it has a face, so a broken body never reads as closed.
    pub fn is_free_edge(&self, edge: EdgeHandle) -> bool {
        self.edge_use(edge) == EdgeUse::Free
    }

    fn edge_use(&self, edge: EdgeHandle) -> EdgeUse {
        let Some(e) = self.edges.get(edge) else {
            return EdgeUse::Free;
        };
        let face_of = |he| {
            self.half_edges
                .get(he)
                .and_then(|h| self.loops.get(h.loop_))
                .map(|l| l.face)
        };
        let Some(he) = self.half_edges.get(e.half_edge) else {
            return EdgeUse::Free;
        };
        match (face_of(e.half_edge), face_of(he.twin)) {
            (None, None) => EdgeUse::Wire,
            (Some(a), Some(b)) if a != b => EdgeUse::Bounded,
            (Some(a), Some(_)) if self.face_closed_across_seam(a) => EdgeUse::Bounded,
            _ => EdgeUse::Free,
        }
    }

    fn edge_shell(&self, edge: EdgeHandle) -> Option<ShellHandle> {
        let e = self.edges.get(edge)?;
        let h = self.half_edges.get(e.half_edge)?;
        let l = self.loops.get(h.loop_)?;
        self.faces.get(l.face).map(|f| f.shell)
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
