//! Honest per-face UV trim windows for the self-intersection pair
//! machinery (PK-328 root fix).
//!
//! The pair scan derived its windows from OUTER-LOOP VERTEX POSITIONS
//! (`exact_bounded_parameter_range`). A closed analytic face carries
//! its loop vertices only at seam anchors — a bore-skin cylinder's
//! outer loop is two seam-anchored circles, so every vertex projects
//! to ONE angle and the u-window collapses to a point. The degenerate-
//! range guard then skipped the pair silently, and a hollow whose
//! void skins genuinely crossed (bore skin r = 0.356 through the side
//! skin at y = 0.156, 0.012 deep) passed `validate_all` clean — the
//! wrong-shape-served-clean class the validator exists to catch.
//! Special arms (cylinder×cylinder, coaxial torus) masked the class
//! for their pairs; the generic arm had it for every closed face.
//!
//! The honest window is the PCURVE-sampled trim window: pcurves are
//! chart-coherent per loop, so a wrapped loop's samples span the full
//! period. Faces without pcurves fall back to the positional window.

use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection;
use crate::topology::{BRep, FaceHandle};

/// Samples per pcurve span when accumulating the window.
const PCURVE_SAMPLES: usize = 8;

/// The face's trim window in UV: the bounding box of its loops'
/// pcurve samples (all loops), falling back to the vertex-positional
/// window when no pcurve contributed.
pub(crate) fn face_trim_uv_window(
    brep: &BRep,
    fh: FaceHandle,
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    let mut lo = (f64::MAX, f64::MAX);
    let mut hi = (f64::MIN, f64::MIN);
    let mut any = false;
    if let Some(face) = brep.faces.get(fh) {
        let loops = std::iter::once(face.outer_loop).chain(face.inner_loops.iter().copied());
        for lh in loops {
            let Some(loop_data) = brep.loops.get(lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            if !brep.half_edges.contains_key(start) {
                continue;
            }
            let mut cursor = start;
            let bound = brep.half_edges.len();
            let mut steps = 0usize;
            while let Some(he) = brep.half_edges.get(cursor) {
                if let Some(pcurve) = he.pcurve.as_ref() {
                    let (t0, t1) = pcurve.parameter_range();
                    for k in 0..=PCURVE_SAMPLES {
                        let t = t0 + (t1 - t0) * k as f64 / PCURVE_SAMPLES as f64;
                        let q = pcurve.point_at(t);
                        lo = (lo.0.min(q.x), lo.1.min(q.y));
                        hi = (hi.0.max(q.x), hi.1.max(q.y));
                        any = true;
                    }
                }
                cursor = he.next;
                steps += 1;
                if cursor == start || steps > bound {
                    break;
                }
            }
        }
    }
    if any {
        ((lo.0, hi.0), (lo.1, hi.1))
    } else {
        let positions = brep.face_outer_loop_positions(fh);
        uv_projection::exact_bounded_parameter_range(&positions, surface)
    }
}
