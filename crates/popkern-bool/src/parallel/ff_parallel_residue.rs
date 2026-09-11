//! Trim-residue record intake for the parallel pair loop (split from
//! `ff_parallel.rs` under the 500-line structural budget).

use popkern_core::topology::FaceHandle;

use crate::intersect::IntersectionRecord;

/// Convert the clip's trim residues into single-face
/// `IntersectionRecord`s — see `IntersectionRecord::residue_owner_a`
/// (PK-465). Residue ends are NOT junction-refined: refining either
/// end re-strands the m6 mouth chamfer (measured, PK-472: 4 sides
/// both ways); the junction-cluster geometry is reconciled
/// post-assembly instead.
pub(super) fn push_residue_records(
    records: &mut Vec<IntersectionRecord>,
    residues: Vec<crate::clip::TrimResidue>,
    face_a: FaceHandle,
    face_b: FaceHandle,
    clip_trace: bool,
) {
    for (owner_is_a, piece, _shared_ends) in residues {
        if clip_trace {
            use popkern_core::geometry::curve::CurveEval as _;
            let (t0, t1) = piece.curve_3d.parameter_range();
            let p0 = piece.curve_3d.point_at(t0);
            let p1 = piece.curve_3d.point_at(t1);
            eprintln!(
                "  [residue] A{face_a:?} x B{face_b:?} owner={} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                if owner_is_a { "A" } else { "B" },
                p0.x,
                p0.y,
                p0.z,
                p1.x,
                p1.y,
                p1.z
            );
        }
        records.push(IntersectionRecord {
            residue_owner_a: Some(owner_is_a),
            face_a,
            face_b,
            curve: piece,
            is_tangent: false,
            is_coincident: false,
            start_boundary: Default::default(),
            end_boundary: Default::default(),
        });
    }
}
