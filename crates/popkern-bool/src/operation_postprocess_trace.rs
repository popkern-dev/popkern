//! Post-processing trace census (split from `operation_postprocess.rs`
//! under the 500-line structural budget): per-pass face/edge/self-twin
//! counts (`boolean` topic) and the between-pass torn-edge census
//! (`uv-scramble` topic, PK-485).

use popkern_core::topology::BRep;

/// Per-stage face/edge census line for the post-processing stages
/// (`boolean` topic) — sliver removal silently deleted a cylinder
/// wall once (PK-150); the census makes any face-count change in the
/// post pipeline visible. Under the `uv-scramble` topic it also runs
/// the torn-edge census between passes (PK-485).
pub(crate) fn trace_face_census(label: &str, brep: &BRep) {
    if popkern_core::trace::trace_enabled("uv-scramble") {
        let tears = brep.find_tessellation_non_manifold_edges(9, 1e-2);
        if !tears.is_empty() {
            eprintln!(
                "  [uv-scramble] {label}: {} torn edge(s), worst {:.3e}",
                tears.len(),
                tears
                    .iter()
                    .map(|r| r.max_divergence)
                    .fold(0.0_f64, f64::max)
            );
        }
    }
    if popkern_core::trace::trace_enabled("boolean") {
        let self_twins = brep
            .half_edges
            .iter()
            .filter(|(heh, he)| he.twin == *heh)
            .count();
        eprintln!(
            "  [post] {label}: {} faces {} edges {} self-twinned",
            brep.faces.len(),
            brep.edges.len(),
            self_twins
        );
        if self_twins > 0 && popkern_core::trace::trace_enabled("spur") {
            for (heh, he) in &brep.half_edges {
                if he.twin != heh {
                    continue;
                }
                let face = brep.loops.get(he.loop_).map(|l| l.face);
                let pos = brep.vertices.get(he.origin).map(|v| v.position);
                eprintln!(
                    "    [post] self-twin {heh:?} face {face:?} loop {:?} at {pos:?}",
                    he.loop_
                );
            }
        }
    }
}
