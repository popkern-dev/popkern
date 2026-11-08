//! Deciding regions from a body's face sets (PK-1102): where a sheet's
//! shell lies once it is known to be open, and the regions of a body
//! written before regions were stored.

use std::collections::BTreeSet;

use slotmap::Key;

use super::{BRep, FaceHandle, Region, RegionHandle, RegionKind, Shell, ShellHandle};

impl BRep {
    /// Move every open solid shell (a shell with a free edge) to the void
    /// it lies in, holding both sides of its faces: a sheet bounds no
    /// volume, so it is no solid region's shell. The void partner its
    /// faces fronted onto goes, and so does its solid region when no
    /// shell is left to it. Called where a producer has finished a body
    /// that may be a sheet (the builder, face removal, sheet booleans).
    ///
    /// A shell is left alone when its faces front onto more than one
    /// shell (faces between regions of a general body), or onto the
    /// peripheral shell of a cavity (a solid whose inner shell is open
    /// is not a sheet; validation reports it). Returns the shells moved.
    pub fn settle_open_shells(&mut self) -> Vec<ShellHandle> {
        let mut settled = Vec::new();
        for shell in self.open_shells() {
            let Some(region) = self.shell_region(shell) else {
                continue;
            };
            if self.region_kind(region) != Some(RegionKind::Solid) {
                continue;
            }
            let fronts: BTreeSet<ShellHandle> = self
                .faces
                .values()
                .filter(|f| f.back == shell && f.front != shell)
                .map(|f| f.front)
                .collect();
            let target = match fronts.iter().next() {
                None => self.outside,
                Some(&partner) if fronts.len() == 1 => {
                    let Some(around) = self.shell_region(partner) else {
                        continue;
                    };
                    let peripheral = self.regions[around].shell == partner;
                    if self.region_kind(around) != Some(RegionKind::Void)
                        || (peripheral && around != self.outside)
                    {
                        continue;
                    }
                    around
                }
                Some(_) => continue,
            };
            for face in self.faces.values_mut() {
                if face.back == shell {
                    face.front = shell;
                }
            }
            self.shells[shell].region = target;
            for partner in fronts {
                if !self
                    .faces
                    .values()
                    .any(|f| f.back == partner || f.front == partner)
                {
                    self.shells.remove(partner);
                }
            }
            if self.regions[region].shell == shell {
                let next = self
                    .shells
                    .iter()
                    .find(|(_, s)| s.region == region)
                    .map_or(ShellHandle::null(), |(h, _)| h);
                self.regions[region].shell = next;
            }
            if !self.shells.values().any(|s| s.region == region) {
                self.regions.remove(region);
            }
            settled.push(shell);
        }
        settled
    }

    /// Make every closed shell lying in a void with both sides of its
    /// faces (a sheet whose last hole a repair closed) the peripheral
    /// shell of a new solid region: its faces keep their backs in it and
    /// front onto a new partner shell in the void it lay in. The
    /// material is behind the faces, so the closed sheet must face out
    /// (validation checks it). The counterpart of
    /// [`BRep::settle_open_shells`], called where an operation's intent
    /// is a solid (hole healing, the boundary-ring caps): a closed
    /// surface separates space, so one shell holding the void on both
    /// sides of it is not a valid cell structure (validation reports it
    /// as a shell of two pieces). Returns the regions made.
    pub fn fill_closed_shells(&mut self) -> Vec<RegionHandle> {
        let open = self.open_shells();
        let candidates: Vec<ShellHandle> = self
            .shells
            .iter()
            .filter(|(h, s)| {
                !open.contains(h)
                    && self.region_kind(s.region) == Some(RegionKind::Void)
                    && self.faces.values().any(|f| f.back == *h)
                    && self
                        .faces
                        .values()
                        .filter(|f| f.back == *h || f.front == *h)
                        .all(|f| f.back == *h && f.front == *h)
            })
            .map(|(h, _)| h)
            .collect();
        let mut made = Vec::new();
        for shell in candidates {
            let around = self.shells[shell].region;
            let anchor = self.shells[shell].face;
            let region = self.regions.insert(Region {
                kind: RegionKind::Solid,
                shell,
            });
            self.shells[shell].region = region;
            let partner = self.shells.insert(Shell::new(anchor, around));
            for face in self.faces.values_mut() {
                if face.back == shell {
                    face.front = partner;
                }
            }
            made.push(region);
        }
        made
    }

    /// Give a body written before regions were stored its regions: the
    /// outside, a partner across each solid shell whose faces have no
    /// front yet (in the outside for a region's outer shell, the
    /// peripheral shell of a new cavity region for an inner one), and
    /// open shells settled as sheets ([`BRep::settle_open_shells`]). A
    /// body that already has its regions is unchanged.
    pub fn derive_regions(&mut self) {
        if !self.regions.contains_key(self.outside) {
            self.outside = self.regions.insert(Region {
                kind: RegionKind::Void,
                shell: ShellHandle::null(),
            });
        }
        let unsided: BTreeSet<ShellHandle> = self
            .faces
            .values()
            .filter(|f| f.front.is_null())
            .map(|f| f.back)
            .filter(|s| self.shells.contains_key(*s))
            .collect();
        for shell in unsided {
            let anchor: FaceHandle = self.shells[shell].face;
            let region = self.shells[shell].region;
            let outer = match self.regions.get_mut(region) {
                // A void's shell holds a sheet, both sides of its faces.
                Some(r) if r.kind == RegionKind::Void => {
                    for face in self.faces.values_mut() {
                        if face.back == shell && face.front.is_null() {
                            face.front = shell;
                        }
                    }
                    continue;
                }
                Some(r) => {
                    if !self.shells.contains_key(r.shell) {
                        r.shell = shell;
                    }
                    r.shell == shell
                }
                // A shell of no live region (an unfinished body) is
                // partnered in the outside.
                _ => true,
            };
            let front = if outer {
                self.shells.insert(Shell::new(anchor, self.outside))
            } else {
                let cavity = self.regions.insert(Region {
                    kind: RegionKind::Void,
                    shell: ShellHandle::null(),
                });
                let front = self.shells.insert(Shell::new(anchor, cavity));
                self.regions[cavity].shell = front;
                front
            };
            for face in self.faces.values_mut() {
                if face.back == shell && face.front.is_null() {
                    face.front = front;
                }
            }
        }
        self.settle_open_shells();
        self.settle_wire_edges();
    }
}
