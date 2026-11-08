//! Region sense (unified topology §5.9): a solid region's faces point
//! out of its material.
//!
//! A face's back side lies in the solid region's shell, so its folded
//! normal (`flag ? n : -n`) must point away from the material: the
//! outer shell then encloses positive folded-normal volume and every
//! cavity shell negative. The volume is the sampled wedge integral of
//! [`BRep::enclosed_folded_volume`]; a shell too thin or open to call is
//! not judged. Faces with both sides in the shell (a sheet inside the
//! solid) bound no volume and are left out.

use crate::error::ValidationError;
use crate::topology::{BRep, FaceHandle};

impl BRep {
    /// Every solid region shell whose folded normals point into the
    /// material ([`ValidationError::RegionSense`]).
    pub fn validate_region_sense_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        for region in self.solid_regions() {
            let outer_shell = self.region_outer_shell(region);
            for shell in self.region_shells(region) {
                let faces: Vec<FaceHandle> = self
                    .shell_faces(shell)
                    .into_iter()
                    .filter(|&f| {
                        self.face_sides(f)
                            .is_some_and(|s| s.back == shell && s.front != shell)
                    })
                    .collect();
                if faces.is_empty() {
                    continue;
                }
                let Some(volume) = self.enclosed_folded_volume(&faces) else {
                    continue;
                };
                let outer = Some(shell) == outer_shell;
                if (volume > 0.0) != outer {
                    errors.push(ValidationError::RegionSense {
                        region,
                        shell,
                        outer,
                        volume,
                    });
                }
            }
        }
        errors
    }
}

#[cfg(test)]
#[path = "region_sense_tests.rs"]
mod tests;
