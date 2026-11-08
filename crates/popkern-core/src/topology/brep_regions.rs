//! Regions and the shells that bound them, read through one place
//! (PK-1102); the edits are in `brep_region_edits`.
//!
//! A region is a connected 3-cell of space, [`RegionKind::Solid`]
//! material or [`RegionKind::Void`] space. Every body has one unbounded
//! void region, [`BRep::outside`]; a cavity is a bounded void region. A
//! shell is a connected set of face sides bounding one region, and a
//! face has two sides, its back (the side a correctly oriented face's
//! normal points away from) and its front, each in a shell
//! ([`Face::back`](super::Face::back), [`Face::front`](super::Face::front)).
//!
//! A closed solid's face has its back in a solid region's shell and its
//! front in the void's shell across it, so the solid's outer boundary is
//! two shells: the solid's and its partner in the outside. A sheet's
//! face has both sides in one shell of the void it lies in. A cavity is
//! a solid region's inner shell, partnered by the peripheral shell of the
//! cavity's void region.
//!
//! The sides record which region lies where; they are topology, not a
//! reading of the normals. Turning a body's faces
//! ([`BRep::reverse_body_orientation`]) turns its geometry only, so a
//! reversed solid is inside out until its faces move into a cavity
//! shell (or turn back).
//!
//! The shells holding some face's back side, [`BRep::back_shells`], are
//! the connected face sets: a cube's one, a hollow ball's two (its
//! outside and its cavity's boundary), a sheet's one. Counts and walks
//! over face sets (Euler–Poincaré, operator count contracts, lumps) read
//! them, never the whole shell arena, which also holds the void sides.

use std::collections::BTreeSet;

use super::{BRep, FaceHandle, RegionHandle, RegionKind, ShellHandle};

/// The two shells a face's sides lie in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct FaceSides {
    /// The shell holding the face's back side.
    pub back: ShellHandle,
    /// The shell holding the face's front side.
    pub front: ShellHandle,
}

impl FaceSides {
    /// Both sides in one shell: a sheet's face.
    pub fn both(shell: ShellHandle) -> Self {
        Self {
            back: shell,
            front: shell,
        }
    }
}

impl BRep {
    /// The shells holding the back sides of faces: the body's connected
    /// face sets (see the module docs), in arena order. A solid region's
    /// shell is one even while it has no face (a body being built, or
    /// one a corruption emptied); a void's shell is one when some face
    /// has its back in it (a sheet's).
    pub fn back_shells(&self) -> Vec<ShellHandle> {
        let backs: BTreeSet<ShellHandle> = self.faces.values().map(|f| f.back).collect();
        self.shells
            .iter()
            .filter(|(h, s)| {
                backs.contains(h) || self.region_kind(s.region) == Some(RegionKind::Solid)
            })
            .map(|(h, _)| h)
            .collect()
    }

    /// How many shells hold some face's back side.
    pub fn back_shell_count(&self) -> usize {
        self.back_shells().len()
    }

    /// The shell holding `face`'s back side.
    pub fn face_back_shell(&self, face: FaceHandle) -> Option<ShellHandle> {
        self.faces
            .get(face)
            .map(|f| f.back)
            .filter(|s| self.shells.contains_key(*s))
    }

    /// The shell holding `face`'s front side.
    pub fn face_front_shell(&self, face: FaceHandle) -> Option<ShellHandle> {
        self.faces
            .get(face)
            .map(|f| f.front)
            .filter(|s| self.shells.contains_key(*s))
    }

    /// The shells `face`'s two sides lie in.
    pub fn face_sides(&self, face: FaceHandle) -> Option<FaceSides> {
        self.faces.get(face).map(|f| FaceSides {
            back: f.back,
            front: f.front,
        })
    }

    /// The faces with a side in `shell`, in arena order: for a solid's
    /// shell the faces with their back in it, for the void's shell
    /// across it the same faces by their fronts, for a sheet's shell its
    /// faces.
    pub fn shell_faces(&self, shell: ShellHandle) -> Vec<FaceHandle> {
        self.faces
            .iter()
            .filter(|(_, f)| f.back == shell || f.front == shell)
            .map(|(h, _)| h)
            .collect()
    }

    /// The region `shell` bounds.
    pub fn shell_region(&self, shell: ShellHandle) -> Option<RegionHandle> {
        self.shells
            .get(shell)
            .map(|s| s.region)
            .filter(|r| self.regions.contains_key(*r))
    }

    /// What fills `region`.
    pub fn region_kind(&self, region: RegionHandle) -> Option<RegionKind> {
        self.regions.get(region).map(|r| r.kind)
    }

    /// The solid regions, in arena order.
    pub fn solid_regions(&self) -> Vec<RegionHandle> {
        self.regions
            .iter()
            .filter(|(_, r)| r.kind == RegionKind::Solid)
            .map(|(h, _)| h)
            .collect()
    }

    /// The bounded void regions (cavities), in arena order.
    pub fn cavity_regions(&self) -> Vec<RegionHandle> {
        self.regions
            .iter()
            .filter(|(h, r)| r.kind == RegionKind::Void && *h != self.outside)
            .map(|(h, _)| h)
            .collect()
    }

    /// A region's outer (peripheral) shell; none for the outside.
    pub fn region_outer_shell(&self, region: RegionHandle) -> Option<ShellHandle> {
        self.regions
            .get(region)
            .map(|r| r.shell)
            .filter(|s| self.shells.contains_key(*s))
    }

    /// Every shell bounding `region`: its outer shell first, then the
    /// rest (a solid's cavities, the outside's shells) in arena order.
    pub fn region_shells(&self, region: RegionHandle) -> Vec<ShellHandle> {
        let outer = self.region_outer_shell(region);
        outer
            .into_iter()
            .chain(
                self.shells
                    .iter()
                    .filter(|(h, s)| s.region == region && Some(*h) != outer)
                    .map(|(h, _)| h),
            )
            .collect()
    }

    /// The face `shell` is anchored at (a face with a side in it), when
    /// that face is live.
    pub fn shell_anchor(&self, shell: ShellHandle) -> Option<FaceHandle> {
        self.shells
            .get(shell)
            .map(|s| s.face)
            .filter(|f| self.faces.contains_key(*f))
    }
}

#[cfg(test)]
#[path = "brep_regions_tests.rs"]
mod tests;
