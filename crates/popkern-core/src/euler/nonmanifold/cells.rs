//! The regions a glued complex of faces encloses (PK-1108, S9b).
//!
//! Faces assembled into one body edge by edge (a STEP shell whose edges
//! some face set uses more than twice, an internal face between two
//! cells) start as sheets in the void: every side of every face in the
//! outside's shells. Once their edges are glued into rings, the sides a
//! wedge connects form the pieces of a cell structure: each closed piece
//! bounds one cell of space. A piece enclosing positive volume (its
//! sides read with the region behind a back side and in front of a front
//! side) is a bounded cell, made a solid region; the rest bound the
//! outside. A face between two cells has its back in one piece and its
//! front in the other, so it bounds both.

use crate::topology::{BRep, FaceHandle, Region, RegionHandle, RegionKind, ShellHandle};

impl BRep {
    /// Give every closed piece of the outside's shells that encloses
    /// positive volume a solid region of its own. The outside's shells
    /// are first separated into their pieces (sides connected through
    /// wedges). Returns the regions made.
    pub fn fill_enclosed_cells(&mut self) -> Vec<RegionHandle> {
        let outside = self.outside;
        let shells: Vec<ShellHandle> = self
            .shells
            .iter()
            .filter(|(_, s)| s.region == outside)
            .map(|(h, _)| h)
            .collect();
        for shell in shells {
            self.separate_shell_pieces(shell);
        }
        let open = self.open_shells();
        let candidates: Vec<ShellHandle> = self
            .shells
            .iter()
            .filter(|(h, s)| s.region == outside && !open.contains(h))
            .map(|(h, _)| h)
            .collect();
        let mut made = Vec::new();
        for shell in candidates {
            let sides: Vec<(FaceHandle, f64)> = self
                .faces
                .iter()
                .filter_map(|(h, f)| match (f.back == shell, f.front == shell) {
                    (true, false) => Some((h, 1.0)),
                    (false, true) => Some((h, -1.0)),
                    _ => None,
                })
                .collect();
            if sides.is_empty() {
                continue;
            }
            if self.enclosed_signed_volume(&sides).is_some_and(|v| v > 0.0) {
                let region = self.regions.insert(Region {
                    kind: RegionKind::Solid,
                    shell,
                });
                self.shells[shell].region = region;
                made.push(region);
            }
        }
        made
    }
}
