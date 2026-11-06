//! Boundary-relative classification of a point against a sheet body (§6.x).
//!
//! Sheet bodies do not enclose a volume, so point-set "inside/outside"
//! classification is undefined. [`side_of_sheet`] classifies a point by
//! which side of the oriented surface it lies on: the "front" side is
//! the half-space the surface normal points into.

use popkern_core::topology::BRep;

/// Side of a sheet body relative to a query point (§6.x boundary-relative).
///
/// Sheet bodies do not enclose a volume, so point-set "inside/outside"
/// classification is undefined. Instead, classify the point by which
/// side of the oriented surface it lies on. The "front" side is the
/// half-space the surface normal points into.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SideOfSheet {
    /// Point is on the side the sheet's outward normal points to.
    Front,
    /// Point is on the opposite side.
    Back,
    /// Point lies on the sheet itself within tolerance.
    On,
    /// Sheet has no orientable face near the query — undefined.
    Indeterminate,
}

/// Classify a query point against a sheet body via boundary-relative
/// semantics.
///
/// Finds the face nearest the point and reports which side of that
/// face's oriented surface the point lies on. A face is near where its
/// surface's closest point lies inside its trims: an unbounded plane's
/// foot past its face's edges is no point of the face (a cube's side
/// wall plane is nearer a point above the cube than its top face is).
/// When no face's foot lies inside its trims (the point is beyond every
/// face's edges), the face with the nearest boundary sample is read at
/// its surface foot. If the point is closer than `tolerance.linear()`,
/// returns `On`.
///
/// This is the building block §6.x asks for: open-shell classification
/// uses "which side of the oriented surface" instead of "inside/outside
/// volume". Callers can use [`side_of_sheet`] in their own
/// sheet-aware boolean wrappers.
pub fn side_of_sheet(
    sheet: &BRep,
    point: &popkern_core::geometry::Point3,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> SideOfSheet {
    use popkern_core::geometry::surface::SurfaceEval;

    // (distance, signed distance) of the nearest face whose foot lies
    // in its trims, and of the face with the nearest boundary sample.
    let mut in_trim: Option<(f64, f64)> = None;
    let mut by_boundary: Option<(f64, f64)> = None;
    for (fh, face) in &sheet.faces {
        let Some(surf) = face.surface.as_ref() else {
            continue;
        };
        let Some((u, v)) = surf.geometry.closest_point(point) else {
            continue;
        };
        let p_surf = surf.geometry.point_at(u, v);
        let mut normal = surf.geometry.normal_at(u, v);
        if !surf.orientation {
            normal = normal * (-1.0);
        }
        let n_len = normal.length();
        if n_len < 1e-30 {
            continue;
        }
        let to_point = *point - p_surf;
        let dist = to_point.length();
        let signed = to_point.dot(&normal) / n_len;
        let foot = popkern_core::geometry::Point3::new(u, v, 0.0);
        if crate::split_inner_loop::point_inside_face_region_uv(sheet, fh, foot, tolerance) {
            if in_trim.is_none_or(|(best, _)| dist < best) {
                in_trim = Some((dist, signed));
            }
            continue;
        }
        let edge_dist = sheet
            .face_boundary_sample_points(fh, 8)
            .iter()
            .map(|q| q.distance_to(point))
            .fold(f64::INFINITY, f64::min);
        if edge_dist.is_finite() && by_boundary.is_none_or(|(best, _)| edge_dist < best) {
            by_boundary = Some((edge_dist, signed));
        }
    }

    // A face reached inside its trims wins unless some face's edge is
    // nearer still (the point then lies beyond that face's edge).
    let best = match (in_trim, by_boundary) {
        (Some(t), Some(b)) if b.0 < t.0 => Some(b),
        (Some(t), _) => Some(t),
        (None, b) => b,
    };
    let Some((best_dist, best_signed_dist)) = best else {
        return SideOfSheet::Indeterminate;
    };

    if best_dist < tolerance.linear() {
        SideOfSheet::On
    } else if best_signed_dist >= 0.0 {
        SideOfSheet::Front
    } else {
        SideOfSheet::Back
    }
}
