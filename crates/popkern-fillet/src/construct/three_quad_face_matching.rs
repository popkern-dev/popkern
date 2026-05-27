//! Sub-face matching helpers for pairing cap patches with fillet boundaries.

use super::*;

pub(super) fn face_has_cylinder_surface(brep: &BRep, face: FaceHandle) -> bool {
    matches!(
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry),
        Some(SurfaceGeometry::Cylinder { .. })
    )
}

pub(super) fn half_edge_endpoints_match_curve(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    curve: &popkern_core::geometry::NurbsCurve,
    tolerance: f64,
) -> bool {
    use popkern_core::geometry::curve::CurveEval;

    let Some(he_data) = brep.half_edges.get(he) else {
        return false;
    };
    let Some(twin_data) = brep.half_edges.get(he_data.twin) else {
        return false;
    };
    let (Some(origin), Some(dest)) = (
        brep.vertices.get(he_data.origin),
        brep.vertices.get(twin_data.origin),
    ) else {
        return false;
    };
    let (t0, t1) = curve.parameter_range();
    let p0 = curve.point_at(t0);
    let p1 = curve.point_at(t1);
    let forward = origin.position.distance_to(&p0) + dest.position.distance_to(&p1);
    let reversed = origin.position.distance_to(&p1) + dest.position.distance_to(&p0);
    forward.min(reversed) <= tolerance
}

pub(super) fn face_surface_matches(
    brep: &BRep,
    face: FaceHandle,
    expected: &SurfaceGeometry,
) -> bool {
    let Some(actual) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    match (&actual.geometry, expected) {
        (
            SurfaceGeometry::Cylinder {
                origin: ao,
                axis: aa,
                radius: ar,
            },
            SurfaceGeometry::Cylinder {
                origin: eo,
                axis: ea,
                radius: er,
            },
        ) => {
            let Some(aa) = aa.normalized() else {
                return false;
            };
            let Some(ea) = ea.normalized() else {
                return false;
            };
            let d = *ao - *eo;
            let along = d.dot(&ea);
            let perp = d - ea * along;
            (ar - er).abs() <= er.abs() * 0.01 + 1.0e-6
                && aa.dot(&ea).abs() >= 0.999
                && perp.length() <= 1.0e-3
        }
        _ => false,
    }
}

/// Return the (origin, axis, radius) of the cylinder that the given
/// face's outer-loop edges share with the most twin faces. Used by
/// the Phase 6c-E pyramid path to disambiguate which per-fillet sphere
/// to assign to each sub-face when all per-fillet radii are equal —
/// matching by axis identity (sphere centre projects onto axis) rather
/// than by radius.
pub(super) fn primary_cylinder_identity_of_face(
    brep: &BRep,
    face: FaceHandle,
) -> Option<(Point3, popkern_core::geometry::Vec3, f64)> {
    let face_data = brep.faces.get(face)?;
    let outer_loop = face_data.outer_loop;
    let loop_data = brep.loops.get(outer_loop)?;
    let start = loop_data.half_edge;

    // Bucket by (origin, axis) fingerprint (rounded). Track the count
    // per bucket; pick the bucket with the highest count.
    type CylKey = (i64, i64, i64, i64, i64, i64);
    type CylBucket = (Point3, popkern_core::geometry::Vec3, f64, usize);
    let mut counts: std::collections::BTreeMap<CylKey, CylBucket> =
        std::collections::BTreeMap::new();
    let mut cursor = start;
    let mut safety = 16usize;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        let he = brep.half_edges.get(cursor)?;
        let twin = he.twin;
        let next = he.next;
        let twin_loop = brep.half_edges.get(twin)?.loop_;
        let twin_face = brep.loops.get(twin_loop)?.face;
        if let Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) = brep
            .faces
            .get(twin_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        {
            let key = (
                (origin.x * 1e6).round() as i64,
                (origin.y * 1e6).round() as i64,
                (origin.z * 1e6).round() as i64,
                (axis.x * 1e6).round() as i64,
                (axis.y * 1e6).round() as i64,
                (axis.z * 1e6).round() as i64,
            );
            counts.entry(key).or_insert((*origin, *axis, *radius, 0)).3 += 1;
        }
        cursor = next;
        if cursor == start {
            break;
        }
    }

    let mut best: Option<(Point3, popkern_core::geometry::Vec3, f64, usize)> = None;
    for &(o, a, r, c) in counts.values() {
        match best {
            Some((_, _, _, bc)) if c <= bc => {}
            _ => best = Some((o, a, r, c)),
        }
    }
    best.map(|(o, a, r, _)| (o, a, r))
}

/// Phase 6c-D Stage 4.7e: build a permutation mapping
/// boundary-segment order (`k = 0..N`) to sub-face index in
/// `sub_faces_raw`.  Returns `Some(perm)` where `sub_faces_raw[perm[k]]`
/// is the sub-face geometrically adjacent to fillet `k`'s
/// cylinder; `None` when matching fails or input shapes are
/// inconsistent.
///
/// Match criterion: the sub-face's PRIMARY cylinder twin (most-
/// frequent cylinder face on its outer loop, via
/// `primary_cylinder_identity_of_face`) shares (origin, axis,
/// radius) with `fillet_surfaces[seg.fillet_index].surface` to
/// numerical tolerance.  Each fillet matches at most one
/// sub-face; sub-faces without a match keep their raw position
/// (best-effort).
pub(super) fn boundary_segment_to_sub_face_permutation(
    brep: &BRep,
    sub_faces_raw: &[FaceHandle],
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<Vec<usize>> {
    let n = sub_faces_raw.len();
    if n == 0 || boundary.segments.len() != 2 * n {
        return None;
    }
    if let Some(perm) = boundary_endpoint_to_sub_face_permutation(brep, sub_faces_raw, boundary) {
        return Some(perm);
    }
    // For each sub-face, identify its primary cyl-twin's params.
    let sub_face_cyls: Vec<Option<(Point3, popkern_core::geometry::Vec3, f64)>> = sub_faces_raw
        .iter()
        .map(|&fh| primary_cylinder_identity_of_face(brep, fh))
        .collect();
    // For each boundary fillet (segment 2k), find the sub-face whose
    // primary cyl-twin matches that fillet's cylinder geometry.
    let mut perm = vec![usize::MAX; n];
    let mut consumed = vec![false; n];
    let mut all_matched = true;
    #[allow(clippy::needless_range_loop)]
    for k in 0..n {
        let seg = &boundary.segments[2 * k];
        let fillet_idx = match seg {
            BoundarySegment::FilletEnd { fillet_index, .. } => *fillet_index,
            _ => return None,
        };
        let fillet_cyl = match &fillet_surfaces.get(fillet_idx)?.surface {
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => (*origin, *axis, *radius),
            _ => {
                // Non-cylinder fillet (Torus etc.) — skip
                // permutation; downstream behaves as before.
                all_matched = false;
                continue;
            }
        };
        let fillet_axis_unit = fillet_cyl.1.normalized()?;
        let mut best: Option<(usize, f64)> = None;
        for (i, sf_cyl_opt) in sub_face_cyls.iter().enumerate() {
            if consumed[i] {
                continue;
            }
            let Some(sf_cyl) = sf_cyl_opt else {
                continue;
            };
            // Radius match.
            if (sf_cyl.2 - fillet_cyl.2).abs() > fillet_cyl.2 * 0.01 + 1e-6 {
                continue;
            }
            // Axis direction match.
            let sf_axis_unit = match sf_cyl.1.normalized() {
                Some(v) => v,
                None => continue,
            };
            let axis_dot = sf_axis_unit.dot(&fillet_axis_unit).abs();
            if axis_dot < 0.999 {
                continue;
            }
            // Origin lies on the fillet's axis line.
            let d = sf_cyl.0 - fillet_cyl.0;
            let along = d.dot(&fillet_axis_unit);
            let perp = d - fillet_axis_unit * along;
            let perp_dist = perp.length();
            if perp_dist > 1e-3 {
                continue;
            }
            // Score: lower is better (more aligned + closer perp).
            let score = perp_dist + (1.0 - axis_dot);
            if best.map(|(_, bs)| score < bs).unwrap_or(true) {
                best = Some((i, score));
            }
        }
        if let Some((i, _)) = best {
            perm[k] = i;
            consumed[i] = true;
        } else {
            all_matched = false;
        }
    }
    if !all_matched {
        return None;
    }
    Some(perm)
}

fn boundary_endpoint_to_sub_face_permutation(
    brep: &BRep,
    sub_faces_raw: &[FaceHandle],
    boundary: &BlendBoundary,
) -> Option<Vec<usize>> {
    let n = sub_faces_raw.len();
    let mut perm = vec![usize::MAX; n];
    let mut consumed = vec![false; n];
    for (k, slot) in perm.iter_mut().enumerate().take(n) {
        let fillet_segment = boundary.segments.get(2 * k)?;
        let face_segment = boundary.segments.get(2 * k + 1)?;
        let mut best: Option<(usize, usize)> = None;
        for (face_idx, &face) in sub_faces_raw.iter().enumerate() {
            if consumed[face_idx] {
                continue;
            }
            let score = usize::from(subface_has_boundary_segment(brep, face, fillet_segment))
                + usize::from(subface_has_boundary_segment(brep, face, face_segment));
            if score == 0 {
                continue;
            }
            if best
                .map(|(_, best_score)| score > best_score)
                .unwrap_or(true)
            {
                best = Some((face_idx, score));
            }
        }
        let (face_idx, _) = best?;
        *slot = face_idx;
        consumed[face_idx] = true;
    }
    Some(perm)
}

fn subface_has_boundary_segment(brep: &BRep, face: FaceHandle, segment: &BoundarySegment) -> bool {
    let curve = match segment {
        BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => curve,
    };
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return false;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..(brep.half_edges.len() + 1) {
        if half_edge_endpoints_match_curve(brep, cursor, curve, 1.0e-5) {
            return true;
        }
        let Some(he) = brep.half_edges.get(cursor) else {
            return false;
        };
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    false
}
