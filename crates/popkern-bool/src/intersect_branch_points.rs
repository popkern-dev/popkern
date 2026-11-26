//! §4.2 branch-point labeling: record endpoints that coincide with
//! another record's endpoint are marked `OnBranchPoint`. Split out of
//! `intersect.rs` to keep the parent under the 500-line structural cap.

use popkern_core::tolerance::Tolerance;

use super::{EndpointBoundaryType, IntersectionRecord};

/// Label endpoints of intersection records that coincide with another
/// curve's endpoint as `OnBranchPoint`.
///
/// Two endpoints "coincide" when their 3D positions are within
/// `tolerance.linear() * 10.0`. The 10× factor matches the consumed
/// proxy used by the SSI marching loop in `popkern-core`, so it
/// catches the same junctions the solver itself treats as coincident.
///
/// The check is O(n²) in the number of intersection records, which is
/// fine for typical boolean operations (well under 100 curves). The
/// label is consumed by future split logic that needs to preserve
/// branch topology when subdividing faces; today it is informational.
pub(super) fn label_branch_point_endpoints(
    records: &mut [IntersectionRecord],
    tolerance: &Tolerance,
) {
    use popkern_core::geometry::CurveEval;
    if records.len() < 2 {
        return;
    }

    let snap = tolerance.linear() * 10.0;
    let snap_sq = snap * snap;

    // Snapshot every endpoint with its (record_idx, is_start) location.
    let endpoints: Vec<(usize, bool, popkern_core::geometry::Point3)> = records
        .iter()
        .enumerate()
        .flat_map(|(i, r)| {
            let (t0, t1) = r.curve.curve_3d.parameter_range();
            let p0 = r.curve.curve_3d.point_at(t0);
            let p1 = r.curve.curve_3d.point_at(t1);
            [(i, true, p0), (i, false, p1)]
        })
        .collect();

    let mut is_branch = vec![[false; 2]; records.len()];

    for i in 0..endpoints.len() {
        for j in (i + 1)..endpoints.len() {
            // Endpoints from the same record are not branches between
            // each other (they are the two ends of one curve).
            if endpoints[i].0 == endpoints[j].0 {
                continue;
            }
            let dx = endpoints[i].2.x - endpoints[j].2.x;
            let dy = endpoints[i].2.y - endpoints[j].2.y;
            let dz = endpoints[i].2.z - endpoints[j].2.z;
            if dx * dx + dy * dy + dz * dz < snap_sq {
                let (idx_i, start_i, _) = endpoints[i];
                let (idx_j, start_j, _) = endpoints[j];
                is_branch[idx_i][if start_i { 0 } else { 1 }] = true;
                is_branch[idx_j][if start_j { 0 } else { 1 }] = true;
            }
        }
    }

    for (i, record) in records.iter_mut().enumerate() {
        if is_branch[i][0] {
            record.start_boundary = EndpointBoundaryType::OnBranchPoint;
        }
        if is_branch[i][1] {
            record.end_boundary = EndpointBoundaryType::OnBranchPoint;
        }
    }
}
