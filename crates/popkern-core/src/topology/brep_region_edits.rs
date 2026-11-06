//! Making and changing regions and shells (PK-1102): the one place
//! outside the construction layer that writes the shell and region
//! arenas. See `brep_regions` for the model.

use std::collections::{BTreeMap, BTreeSet};

use slotmap::Key;

use super::{BRep, FaceHandle, FaceSides, Region, RegionHandle, RegionKind, Shell, ShellHandle};

impl BRep {
    /// Anchor `shell` at `face`.
    pub fn set_shell_anchor(&mut self, shell: ShellHandle, face: FaceHandle) {
        if let Some(s) = self.shells.get_mut(shell) {
            s.face = face;
        }
    }

    /// Anchor both of `sides`' shells at `face`.
    pub fn set_sides_anchor(&mut self, sides: FaceSides, face: FaceHandle) {
        self.set_shell_anchor(sides.back, face);
        self.set_shell_anchor(sides.front, face);
    }

    /// Anchor each shell `face`'s sides lie in at `face` when the shell
    /// has no live anchor yet (the first face made in a new shell).
    pub fn anchor_unset_sides(&mut self, face: FaceHandle) {
        let Some(sides) = self.face_sides(face) else {
            return;
        };
        for shell in [sides.back, sides.front] {
            if self.shell_anchor(shell).is_none() {
                self.set_shell_anchor(shell, face);
            }
        }
    }

    /// Put `face`'s two sides in `sides`' shells.
    pub fn set_face_sides(&mut self, face: FaceHandle, sides: FaceSides) {
        if let Some(f) = self.faces.get_mut(face) {
            f.back = sides.back;
            f.front = sides.front;
        }
    }

    /// The sides a new face added to `shell` takes, read off a face
    /// already there: a solid shell's faces front onto its partner, a
    /// sheet's lie in it on both sides. `None` when no face has its back
    /// in `shell` and its region is solid (an empty solid shell's
    /// partner is known only to whoever made it).
    pub fn sides_beside(&self, shell: ShellHandle) -> Option<FaceSides> {
        if let Some(f) = self.faces.values().find(|f| f.back == shell) {
            return Some(FaceSides {
                back: shell,
                front: f.front,
            });
        }
        match self.region_kind(self.shell_region(shell)?)? {
            RegionKind::Void => Some(FaceSides::both(shell)),
            RegionKind::Solid => None,
        }
    }

    /// Before `face` goes, move the anchor of each shell its sides lie
    /// in off it onto another face of that shell. `false` when such a
    /// shell has no other face (the caller decides: the shell goes too,
    /// or the removal is refused); `true` when every anchor was
    /// elsewhere or has moved.
    pub fn move_shell_anchor_off(&mut self, face: FaceHandle) -> bool {
        let Some(sides) = self.face_sides(face) else {
            return true;
        };
        let mut moved = true;
        for shell in [sides.back, sides.front] {
            if self.shells.get(shell).is_none_or(|s| s.face != face) {
                continue;
            }
            match self
                .shell_faces(shell)
                .into_iter()
                .find(|&candidate| candidate != face)
            {
                Some(replacement) => self.shells[shell].face = replacement,
                None => moved = false,
            }
        }
        moved
    }

    /// Re-anchor every shell whose anchor face is gone onto a live face
    /// with a side in it. Returns the shells left with no face at all
    /// (the caller removes them or declines).
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

    /// Make a solid region bounded by one new shell anchored at `face`
    /// (null while the shell has no face yet), its partner across it a
    /// new shell of the outside. Returns the region and the sides its
    /// faces take.
    pub fn add_solid_region(&mut self, face: FaceHandle) -> (RegionHandle, FaceSides) {
        let region = self.regions.insert(Region {
            kind: RegionKind::Solid,
            shell: ShellHandle::null(),
        });
        let back = self.shells.insert(Shell { face, region });
        let front = self.shells.insert(Shell {
            face,
            region: self.outside,
        });
        self.regions[region].shell = back;
        (region, FaceSides { back, front })
    }

    /// Add a shell to `region`, anchored at `face`, and return the sides
    /// its faces take. A solid region's first shell is its outer one,
    /// partnered in the outside; a later one bounds a cavity, partnered
    /// by the peripheral shell of a new void region. A void region's
    /// shell holds a sheet, both sides of its faces.
    pub fn add_region_shell(&mut self, region: RegionHandle, face: FaceHandle) -> FaceSides {
        let back = self.shells.insert(Shell { face, region });
        if self.region_kind(region) != Some(RegionKind::Solid) {
            return FaceSides::both(back);
        }
        let front = if self.region_outer_shell(region).is_none() {
            self.regions[region].shell = back;
            self.shells.insert(Shell {
                face,
                region: self.outside,
            })
        } else {
            let cavity = self.regions.insert(Region {
                kind: RegionKind::Void,
                shell: ShellHandle::null(),
            });
            let front = self.shells.insert(Shell {
                face,
                region: cavity,
            });
            self.regions[cavity].shell = front;
            front
        };
        FaceSides { back, front }
    }

    /// Turn `regions` inside out into `host`: the regions of a body
    /// whose faces were reversed ([`BRep::reverse_body_orientation`])
    /// and which now lies inside `host`, as a subtracted tool enclosed
    /// by the target. Every face with a side in one of their shells
    /// swaps sides, each of them swaps kind (the tool's material
    /// becomes the host's cavity, a cavity of the tool the material
    /// left inside it), and their shells in the outside move to `host`.
    pub fn invert_regions_into(&mut self, regions: &[RegionHandle], host: RegionHandle) {
        let inverted: BTreeSet<RegionHandle> = regions.iter().copied().collect();
        let inner: BTreeSet<ShellHandle> = self
            .shells
            .iter()
            .filter(|(_, s)| inverted.contains(&s.region))
            .map(|(h, _)| h)
            .collect();
        let mut outer: BTreeSet<ShellHandle> = BTreeSet::new();
        for face in self.faces.values_mut() {
            if inner.contains(&face.back) || inner.contains(&face.front) {
                std::mem::swap(&mut face.back, &mut face.front);
                for side in [face.back, face.front] {
                    if !inner.contains(&side) {
                        outer.insert(side);
                    }
                }
            }
        }
        for r in &inverted {
            if let Some(region) = self.regions.get_mut(*r) {
                region.kind = match region.kind {
                    RegionKind::Solid => RegionKind::Void,
                    RegionKind::Void => RegionKind::Solid,
                };
            }
        }
        for shell in outer {
            if self
                .shells
                .get(shell)
                .is_some_and(|s| s.region == self.outside)
            {
                self.shells[shell].region = host;
            }
        }
    }

    /// Keep only the solid regions in `keep`. Every other solid region
    /// goes with its shells and the faces whose back lies in them; the
    /// void shells those faces fronted onto go when no face is left in
    /// them, and a cavity of a dropped region goes with it, any body
    /// floating in the cavity moving to the void around the dropped
    /// region. Loops, fins, edges and vertices are the caller's (a
    /// measurement view needs only faces; piece extraction removes the
    /// rest by what the kept faces use).
    pub fn retain_solid_regions(&mut self, keep: &[RegionHandle]) {
        let dropped: Vec<RegionHandle> = self
            .solid_regions()
            .into_iter()
            .filter(|r| !keep.contains(r))
            .collect();
        // The void each dropped region sits in, and the cavities it
        // encloses, read before anything goes.
        let mut moved_to: BTreeMap<RegionHandle, RegionHandle> = BTreeMap::new();
        for &r in &dropped {
            let around = self.void_across(self.region_outer_shell(r));
            for inner in self.region_shells(r).into_iter().skip(1) {
                let cavity = self.void_across(Some(inner));
                if cavity != self.outside && cavity != around {
                    moved_to.insert(cavity, around);
                }
            }
        }
        let resolve = |mut r: RegionHandle| {
            while let Some(&next) = moved_to.get(&r) {
                r = next;
            }
            r
        };
        let shells: BTreeSet<ShellHandle> = dropped
            .iter()
            .flat_map(|&r| self.region_shells(r))
            .collect();
        let sided = |b: &BRep| -> BTreeSet<ShellHandle> {
            b.faces.values().flat_map(|f| [f.back, f.front]).collect()
        };
        let sided_before = sided(self);
        self.faces.retain(|_, f| !shells.contains(&f.back));
        let sided_after = sided(self);
        let cavities: BTreeSet<RegionHandle> = moved_to.keys().copied().collect();
        for shell in self.shells.values_mut() {
            if cavities.contains(&shell.region) {
                shell.region = resolve(shell.region);
            }
        }
        // A shell goes with its region, or when its last face went here.
        self.shells.retain(|h, _| {
            !shells.contains(&h) && (sided_after.contains(&h) || !sided_before.contains(&h))
        });
        self.regions
            .retain(|h, _| !dropped.contains(&h) && !cavities.contains(&h));
    }

    /// The void region across `shell` (the region its faces' fronts lie
    /// in); the outside when it has no faces or they front onto no live
    /// shell.
    fn void_across(&self, shell: Option<ShellHandle>) -> RegionHandle {
        shell
            .and_then(|s| self.faces.values().find(|f| f.back == s))
            .and_then(|f| self.shell_region(f.front))
            .unwrap_or(self.outside)
    }
}
