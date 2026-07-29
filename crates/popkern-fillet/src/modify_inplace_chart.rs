//! Strip-chart pcurve re-authoring for in-place blend modification.
//!
//! An in-place radius modify swaps the strip face's surface, which
//! retires the old chart: every strip-side pcurve on the modified
//! face is then stale in a chart that no longer exists. A
//! constant->constant swap survives only because both cylinder
//! charts share the angle/station convention; a variable sweep
//! chart does not — PK-179 measured terminal settlement sampling
//! the new sweep at old cylinder angles, landing 2.1 off the body
//! and declining every corner. So the face's half-edge pcurves and
//! face-use overrides are cleared and recomputed against the NEW
//! surface before any settlement reads them.

use popkern_core::topology::{BRep, FaceHandle, FaceUseTrimKey};

use crate::error::BlendError;

/// Clear and recompute every pcurve on `strip_face` against its
/// (just swapped) surface. Transactional in effect: recompute
/// failure surfaces as a structured `BlendError` and the caller
/// declines the modify.
pub(crate) fn reauthor_strip_chart_pcurves(
    out: &mut BRep,
    strip_face: FaceHandle,
) -> Result<(), BlendError> {
    let strip_half_edges: Vec<_> = {
        let Some(face) = out.faces.get(strip_face) else {
            return Ok(());
        };
        let mut loops = vec![face.outer_loop];
        loops.extend(face.inner_loops.iter().copied());
        let mut hes = Vec::new();
        for lh in loops {
            let Some(loop_data) = out.loops.get(lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut cursor = start;
            loop {
                hes.push(cursor);
                let Some(he) = out.half_edges.get(cursor) else {
                    break;
                };
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
        }
        hes
    };
    for he in strip_half_edges {
        if let Some(he_data) = out.half_edges.get_mut(he) {
            he_data.pcurve = None;
        }
        out.face_use_pcurves.remove(&FaceUseTrimKey {
            face: strip_face,
            half_edge: he,
        });
    }
    out.recompute_face_pcurves(strip_face)
        .map_err(|e| BlendError::InvalidParameter {
            parameter: "strip_face.pcurves".to_string(),
            value: format!("{strip_face:?}"),
            reason: format!("re-authoring strip pcurves in the new chart failed: {e}"),
        })
}
