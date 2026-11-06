//! Regions and the shells that bound them, read and made through one
//! place (PK-1102).
//!
//! A region is a connected 3-cell of space: a body's solid material or
//! a void (the unbounded outside, a cavity). A shell is a connected set
//! of face sides bounding one region, and a face has two sides, its
//! back (the side its effective normal points away from) and its front.
//! A closed solid's faces have their backs in the solid region's
//! shells; a sheet's faces lie in one shell with both their sides.
//!
//! The shells holding some face's back side, [`BRep::back_shells`], are
//! the connected face sets: a cube's one, a hollow ball's two (its
//! outside and its cavity's boundary), a sheet's one. Counts and walks
//! over face sets (Euler–Poincaré, operator count contracts, lumps) read
//! them, never the whole shell arena, which also holds the void sides.

use slotmap::Key;

use super::{BRep, FaceHandle, Shell, ShellHandle, Solid, SolidHandle};

impl BRep {
    /// The shells holding some face's back side: the body's connected
    /// face sets (see the module docs), in arena order.
    pub fn back_shells(&self) -> Vec<ShellHandle> {
        self.shells.keys().collect()
    }

    /// How many shells hold some face's back side.
    pub fn back_shell_count(&self) -> usize {
        self.shells.len()
    }

    /// The shell holding `face`'s back side.
    pub fn face_back_shell(&self, face: FaceHandle) -> Option<ShellHandle> {
        self.faces
            .get(face)
            .map(|f| f.shell)
            .filter(|s| self.shells.contains_key(*s))
    }

    /// The faces whose back side lies in `shell`, in arena order.
    pub fn shell_faces(&self, shell: ShellHandle) -> Vec<FaceHandle> {
        self.faces
            .iter()
            .filter(|(_, f)| f.shell == shell)
            .map(|(h, _)| h)
            .collect()
    }

    /// The region `shell` bounds.
    pub fn shell_region(&self, shell: ShellHandle) -> Option<SolidHandle> {
        self.shells
            .get(shell)
            .map(|s| s.solid)
            .filter(|r| self.solids.contains_key(*r))
    }

    /// The solid regions, in arena order.
    pub fn solid_regions(&self) -> Vec<SolidHandle> {
        self.solids.keys().collect()
    }

    /// A region's outer (peripheral) shell.
    pub fn region_outer_shell(&self, region: SolidHandle) -> Option<ShellHandle> {
        self.solids
            .get(region)
            .map(|r| r.shell)
            .filter(|s| self.shells.contains_key(*s))
    }

    /// Every shell bounding `region`: its outer shell first, then the
    /// rest (cavities) in arena order.
    pub fn region_shells(&self, region: SolidHandle) -> Vec<ShellHandle> {
        let outer = self.region_outer_shell(region);
        outer
            .into_iter()
            .chain(
                self.shells
                    .iter()
                    .filter(|(h, s)| s.solid == region && Some(*h) != outer)
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

    /// Anchor `shell` at `face`.
    pub fn set_shell_anchor(&mut self, shell: ShellHandle, face: FaceHandle) {
        if let Some(s) = self.shells.get_mut(shell) {
            s.face = face;
        }
    }

    /// Before `face` goes, move its back shell's anchor off it onto
    /// another face of that shell. `false` when the shell has no other
    /// face (the caller decides: the shell goes too, or the removal is
    /// refused); `true` when the anchor was elsewhere or has moved.
    pub fn move_shell_anchor_off(&mut self, face: FaceHandle) -> bool {
        let Some(shell) = self.face_back_shell(face) else {
            return true;
        };
        if self.shells[shell].face != face {
            return true;
        }
        match self
            .faces
            .iter()
            .find(|(candidate, data)| *candidate != face && data.shell == shell)
        {
            Some((replacement, _)) => {
                self.shells[shell].face = replacement;
                true
            }
            None => false,
        }
    }

    /// Re-anchor every shell whose anchor face is gone onto a live face
    /// of the shell. Returns the shells left with no face at all (the
    /// caller removes them or declines).
    pub fn repair_shell_anchors(&mut self) -> Vec<ShellHandle> {
        let mut empty = Vec::new();
        for shell in self.shells.keys().collect::<Vec<_>>() {
            if self.shell_anchor(shell).is_some() {
                continue;
            }
            match self.shell_faces(shell).first() {
                Some(&face) => self.shells[shell].face = face,
                None => empty.push(shell),
            }
        }
        empty
    }

    /// Move `shell` to bound `region` instead (a body's shell becoming
    /// a cavity of the region enclosing it). The region it leaves keeps
    /// its other shells; its outer shell moves to one of them, or none.
    pub fn move_shell_to_region(&mut self, shell: ShellHandle, region: SolidHandle) {
        let Some(from) = self.shell_region(shell) else {
            return;
        };
        if from == region || !self.solids.contains_key(region) {
            return;
        }
        self.shells[shell].solid = region;
        if self.solids[from].shell == shell {
            let next = self
                .shells
                .iter()
                .find(|(_, s)| s.solid == from)
                .map_or(ShellHandle::null(), |(h, _)| h);
            self.solids[from].shell = next;
        }
        if !self.shells.contains_key(self.solids[region].shell) {
            self.solids[region].shell = shell;
        }
    }

    /// Remove `region` when no shell bounds it any longer. Returns
    /// whether it went.
    pub fn remove_empty_region(&mut self, region: SolidHandle) -> bool {
        if self.shells.values().any(|s| s.solid == region) {
            return false;
        }
        self.solids.remove(region).is_some()
    }

    /// Keep only the solid regions in `keep`: every other solid region
    /// goes with the shells bounding it and the faces whose back side
    /// lies in them. Loops, fins, edges and vertices are the caller's
    /// (a measurement view needs only faces; piece extraction removes
    /// the rest by what the kept faces use).
    pub fn retain_solid_regions(&mut self, keep: &[SolidHandle]) {
        let dropped: Vec<SolidHandle> = self
            .solid_regions()
            .into_iter()
            .filter(|r| !keep.contains(r))
            .collect();
        let shells: Vec<ShellHandle> = dropped
            .iter()
            .flat_map(|&r| self.region_shells(r))
            .collect();
        self.faces.retain(|_, f| !shells.contains(&f.shell));
        self.shells.retain(|h, _| !shells.contains(&h));
        self.solids.retain(|h, _| !dropped.contains(&h));
    }

    /// Make a solid region bounded by one new shell anchored at `face`
    /// (null while the shell has no face yet). Returns the region and
    /// its shell.
    pub fn add_solid_region(&mut self, face: FaceHandle) -> (SolidHandle, ShellHandle) {
        let region = self.solids.insert(Solid {
            shell: ShellHandle::null(),
        });
        let shell = self.shells.insert(Shell {
            face,
            solid: region,
        });
        self.solids[region].shell = shell;
        (region, shell)
    }

    /// Add a further shell to `region` (a cavity's boundary), anchored
    /// at `face`. The region's outer shell is unchanged unless it had
    /// none.
    pub fn add_region_shell(&mut self, region: SolidHandle, face: FaceHandle) -> ShellHandle {
        let shell = self.shells.insert(Shell {
            face,
            solid: region,
        });
        if let Some(r) = self.solids.get_mut(region)
            && !self.shells.contains_key(r.shell)
        {
            r.shell = shell;
        }
        shell
    }
}

#[cfg(test)]
#[path = "brep_regions_tests.rs"]
mod tests;
