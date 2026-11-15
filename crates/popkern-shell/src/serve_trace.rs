//! The serve gate's `shell-serve` trace, split from `lib.rs` (500-line
//! cap).

use popkern_core::topology::BRep;

/// The `shell-serve` trace of a result the serve gate rejects: its
/// errors and every edge's curve and ends (the gate's own diagnostic).
pub(crate) fn serve_trace(
    result: &BRep,
    validation_errors: &[popkern_core::error::ValidationError],
) -> Result<(), crate::ShellError> {
    if popkern_core::trace::trace_enabled("shell-serve") {
        for e in validation_errors {
            eprintln!("[serve-trace] {e}");
        }
        for (eh, edge) in &result.edges {
            let ends = {
                let he = &result.half_edges[edge.half_edge];
                let tw = &result.half_edges[he.require_mate()?];
                (
                    result.vertices[he.origin].position,
                    result.vertices[tw.origin].position,
                )
            };
            eprintln!(
                "[serve-trace] {eh:?} curve={:?} ends=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                edge.curve.as_ref().map(|c| format!("{:?}", c.geometry)
                    .chars()
                    .take(90)
                    .collect::<String>()),
                ends.0.x,
                ends.0.y,
                ends.0.z,
                ends.1.x,
                ends.1.y,
                ends.1.z
            );
        }
    }
    Ok(())
}
