//! Inner-loop winding: a face's holes wind against its outer loop.
//!
//! Every loop keeps the face's material on its left (walking the
//! half-edges with the flag-folded normal up), so with the flag `true`
//! the outer loop runs counter-clockwise in UV and each inner loop
//! clockwise. The winding-honesty audit measures the outer loop only.
//! An inner loop wound like the outer one still twins with its
//! neighbours, so the shell closes and every structural check passes,
//! but the material lies on opposite sides of that loop (measured: a
//! mixed-sense scoped union served half inside out through a boss
//! bottom's hole, signed volume -606.7 for +617.9; PK-916, PK-923).
//!
//! The check measures each inner loop with the audit's own shoelace
//! (`loop_signed_uv_area`): loops that wrap a period or touch a pole
//! have no UV winding there and are skipped, as are loops under the
//! audit's perimeter-relative degeneracy floor.

use crate::error::ValidationError;
use crate::topology::BRep;

impl BRep {
    /// One [`ValidationError::InnerLoopWinding`] per determinate inner
    /// loop that winds the same way as its face's flag says the outer
    /// loop does.
    pub fn validate_inner_loop_winding_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        for (face, fd) in &self.faces {
            let Some(sof) = fd.surface.as_ref() else {
                continue;
            };
            for &l in &fd.inner_loops {
                let Ok(area) = self.loop_signed_uv_area(l, &sof.geometry) else {
                    continue;
                };
                if area.signed.abs() <= area.perimeter.powi(2) * 1.0e-12 {
                    continue;
                }
                if (area.signed > 0.0) == sof.orientation {
                    errors.push(ValidationError::InnerLoopWinding {
                        face,
                        loop_: l,
                        signed_uv_area: area.signed,
                        orientation: sof.orientation,
                    });
                }
            }
        }
        errors
    }
}
