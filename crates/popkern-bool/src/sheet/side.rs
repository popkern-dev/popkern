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
/// Finds the closest face on the sheet, evaluates its surface and
/// outward normal at the closest UV, and reports which side the
/// point lies on. If the point is closer than `tolerance.linear()`,
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

    let mut best_dist = f64::INFINITY;
    let mut best_signed_dist = 0.0_f64;
    let mut best_found = false;

    for (_fh, face) in &sheet.faces {
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
        if dist < best_dist {
            best_dist = dist;
            let signed = to_point.dot(&normal) / n_len;
            best_signed_dist = signed;
            best_found = true;
        }
    }

    if !best_found {
        return SideOfSheet::Indeterminate;
    }

    if best_dist < tolerance.linear() {
        SideOfSheet::On
    } else if best_signed_dist >= 0.0 {
        SideOfSheet::Front
    } else {
        SideOfSheet::Back
    }
}
