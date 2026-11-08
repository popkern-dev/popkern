//! make_face_region / kill_face_region: a closed surface separates
//! space, and the space on one side of it becomes a region of its own;
//! and back.
//!
//! A shell holding both sides of every face of a closed surface (a
//! closed sheet assembled face by face, a solid whose region was
//! killed) bounds nothing yet: validation reads it as two pieces, one
//! per side. Making the face region gives one side to a new region of
//! the kind asked (a solid filling a closed sheet, a void bubble inside
//! a solid, the void inside a closed sheet); the other side stays in
//! the region the surface lay in, in a shell of its own. Which side is
//! the caller's: the faces' normals say which way they face, and
//! validation checks a solid region's faces point out of it. Killing a
//! face region gives its side of every face back to the shell across.

use super::EulerError;
use crate::topology::{BRep, FaceSide, Region, RegionHandle, RegionKind, Shell, ShellHandle};

impl BRep {
    /// Give `side` of every face of the closed surface `shell` holds to
    /// a new region of `kind`, the other side staying in a new shell of
    /// the shell's region. Returns the region. Counts: shells +1,
    /// regions +1.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotAClosedSurface`] when the shell holds a face
    /// with only one side in it, a free edge, wire edges or an acorn.
    pub fn make_face_region(
        &mut self,
        shell: ShellHandle,
        side: FaceSide,
        kind: RegionKind,
    ) -> Result<RegionHandle, EulerError> {
        let bad = EulerError::NotAClosedSurface { shell };
        let data = self.shells.get(shell).ok_or(bad.clone())?;
        if !data.wire_edges.is_empty() || data.acorn.is_some() {
            return Err(bad);
        }
        let faces = self.shell_faces(shell);
        let closed_lamina = !faces.is_empty()
            && faces
                .iter()
                .all(|&f| self.faces[f].back == shell && self.faces[f].front == shell);
        if !closed_lamina || self.open_shells().contains(&shell) {
            return Err(bad);
        }
        let around = data.region;
        let anchor = data.face;
        let region = self.regions.insert(Region { kind, shell });
        self.shells[shell].region = region;
        let partner = self.shells.insert(Shell::new(anchor, around));
        for f in faces {
            let face = &mut self.faces[f];
            match side {
                FaceSide::Back => face.front = partner,
                FaceSide::Front => face.back = partner,
            }
        }
        Ok(region)
    }

    /// Give `region`'s side of every face it is bounded by back to the
    /// shell across, the inverse of [`BRep::make_face_region`]: the
    /// region and its shell go, and the surface lies with both sides in
    /// that shell. Counts: shells −1, regions −1.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotAFaceRegion`] for the outside, a region with
    /// other than one shell or one holding wires or an acorn, or one
    /// whose faces' other sides do not all lie in one shell of another
    /// region.
    pub fn kill_face_region(&mut self, region: RegionHandle) -> Result<(), EulerError> {
        let bad = EulerError::NotAFaceRegion { region };
        if region == self.outside || !self.regions.contains_key(region) {
            return Err(bad);
        }
        let shells = self.region_shells(region);
        let [shell] = shells[..] else {
            return Err(bad);
        };
        if !self.shells[shell].wire_edges.is_empty() || self.shells[shell].acorn.is_some() {
            return Err(bad);
        }
        let faces = self.shell_faces(shell);
        let mut across = None;
        for &f in &faces {
            let face = &self.faces[f];
            let other = match (face.back == shell, face.front == shell) {
                (true, false) => face.front,
                (false, true) => face.back,
                _ => return Err(bad),
            };
            if *across.get_or_insert(other) != other {
                return Err(bad);
            }
        }
        let Some(across) = across else {
            return Err(bad);
        };
        if self.shells.get(across).is_none_or(|s| s.region == region) {
            return Err(bad);
        }
        for f in faces {
            let face = &mut self.faces[f];
            if face.back == shell {
                face.back = across;
            }
            if face.front == shell {
                face.front = across;
            }
        }
        self.shells.remove(shell);
        self.regions.remove(region);
        Ok(())
    }
}
