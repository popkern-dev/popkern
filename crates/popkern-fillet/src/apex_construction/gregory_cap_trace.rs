//! Trace helpers for Gregory cap boundary assembly.

use popkern_core::geometry::CurveEval;

use crate::blend::boundary::BoundarySegment;

pub(super) fn trace_gregory_boundary_segments(segments: &[BoundarySegment]) {
    if !crate::direct_apex_trace::enabled() {
        return;
    }
    for (index, segment) in segments.iter().enumerate() {
        let (label, curve) = match segment {
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => (format!("FilletEnd[{fillet_index}]"), curve),
            BoundarySegment::FaceEdge { face, curve } => (format!("FaceEdge[{face:?}]"), curve),
        };
        let (t0, t1) = curve.parameter_range();
        let start = curve.point_at(t0);
        let end = curve.point_at(t1);
        eprintln!(
            "[direct-apex-trace §DirectApex.4/gregory-boundary] segment[{index}] {label}: \
             ({:.6},{:.6},{:.6}) -> ({:.6},{:.6},{:.6})",
            start.x, start.y, start.z, end.x, end.y, end.z
        );
    }
}
