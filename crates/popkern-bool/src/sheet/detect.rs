//! Sheet-body detection and geometric extraction helpers.
//!
//! [`is_sheet_body`] is the entry point used by the boolean engine's
//! pre-check to reject open shells upfront. [`is_planar_sheet`] is the
//! fast-path detector used by the multi-piece dispatcher to identify
//! planar sheets that can be split via half-space boolean operations.
//! [`is_approximately_planar_sheet`] extends that fast path to sheets
//! with non-planar surface types (NURBS, sphere caps) whose vertices
//! happen to lie within a tight envelope of a common best-fit plane.

use popkern_core::topology::BRep;

/// Returns `true` if `brep` is a sheet body (open shell).
///
/// Detection: an edge whose two half-edges share the same face is
/// "free" — but only if the underlying surface is NOT periodic in
/// the direction of that edge. Periodic surfaces (cylinder seam,
/// sphere meridian, torus seams) legitimately have face_a == face_b
/// at their seam edges, and the body is still closed.
///
/// The check therefore: walks every edge, identifies free edges by
/// the same-face criterion, and skips those whose hosting surface
/// reports periodicity. Anything left over is a real open boundary.
///
/// Empty BReps (no edges, no faces) are NOT considered sheet bodies
/// — they are simply degenerate, and the existing pre-check rejects
/// them via the "no faces" finding.
pub fn is_sheet_body(brep: &BRep) -> bool {
    use popkern_core::geometry::surface::SurfaceGeometry;

    if brep.edges.is_empty() {
        return false;
    }

    for (_, edge) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            return true;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            return true;
        };
        let face_a = brep.loops.get(he.loop_).map(|l| l.face);
        let face_b = brep.loops.get(twin.loop_).map(|l| l.face);
        let (face_a, face_b) = match (face_a, face_b) {
            (Some(fa), Some(fb)) => (fa, fb),
            _ => return true,
        };

        if face_a != face_b {
            continue;
        }

        let Some(face) = brep.faces.get(face_a) else {
            return true;
        };
        let Some(surf_on_face) = face.surface.as_ref() else {
            return true;
        };
        let is_periodic = match &surf_on_face.geometry {
            SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Cone { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Torus { .. } => true,
            // A seam edge on a NURBS surface is legitimate when the
            // surface is geometrically closed in one direction
            // (e.g. a swept circular profile — the sweep emits a
            // closed NURBS side face with a seam, not an analytic
            // cylinder).
            SurfaceGeometry::Nurbs(_) => {
                use popkern_core::geometry::surface::SurfaceEval;
                let sfc = &surf_on_face.geometry;
                let ((u0, u1), (v0, v1)) = sfc.parameter_range();
                let gate = brep.tolerance.linear().max(1.0e-9) * 10.0;
                let closed_u = (0..=4).all(|i| {
                    let v = v0 + (v1 - v0) * i as f64 / 4.0;
                    sfc.point_at(u0, v).distance_to(&sfc.point_at(u1, v)) < gate
                });
                let closed_v = (0..=4).all(|i| {
                    let u = u0 + (u1 - u0) * i as f64 / 4.0;
                    sfc.point_at(u, v0).distance_to(&sfc.point_at(u, v1)) < gate
                });
                closed_u || closed_v
            }
            _ => false,
        };
        if !is_periodic {
            return true;
        }
    }
    false
}

/// Returns `true` if `brep` is a closed solid (no free edges).
///
/// Inverse of [`is_sheet_body`] for non-empty BReps. An empty BRep
/// (no edges) is reported as closed because there is nothing to
/// breach the closure condition.
pub fn is_closed_solid(brep: &BRep) -> bool {
    !is_sheet_body(brep)
}

/// Detect whether `sheet` is "effectively planar" and extract its plane.
///
/// Returns `Some((origin, outward_normal))` when every face of the
/// sheet lies on the same supporting plane within `tolerance.linear()`
/// — either because the sheet has a single planar face, or because it
/// has multiple planar faces that all share the same normal and offset
/// along that normal. The returned normal is the face's outward-facing
/// normal (orientation flag applied). Returns `None` for sheets with
/// non-planar surfaces, mixed normals, missing surface geometry, or
/// zero faces.
///
/// Used by [`crate::sheet::boolean_with_sheets`] to dispatch the
/// planar fast path for `SolidMinusSheet` — the resulting plane is
/// passed to [`crate::sheet::split_solid_by_plane`] to produce
/// multiple closed sub-solids via half-space split.
pub fn is_planar_sheet(
    sheet: &BRep,
) -> Option<(popkern_core::geometry::Point3, popkern_core::geometry::Vec3)> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    if sheet.faces.is_empty() {
        return None;
    }

    // Tolerance floor for "same plane" checks. We compare the signed
    // offset of each face's own plane origin along the reference
    // normal — if every face agrees within ~1e-6, treat them as a
    // single plane.
    let offset_tol = sheet.tolerance.linear();
    // Cosine threshold for "same normal" (normals within ~1e-4 rad).
    let angle_tol = 1.0e-4;

    let mut ref_origin: Option<popkern_core::geometry::Point3> = None;
    let mut ref_normal: Option<popkern_core::geometry::Vec3> = None;

    for face in sheet.faces.values() {
        let surf_on_face = face.surface.as_ref()?;
        let SurfaceGeometry::Plane { origin, normal, .. } = &surf_on_face.geometry else {
            return None;
        };
        let face_normal = if surf_on_face.orientation {
            *normal
        } else {
            *normal * -1.0
        };

        match (ref_origin, ref_normal) {
            (None, None) => {
                ref_origin = Some(*origin);
                ref_normal = Some(face_normal);
            }
            (Some(ref_origin_pt), Some(rn)) => {
                // Compare normals (absolute value — opposite-facing
                // planar faces are still on the same plane; we keep
                // the first-seen direction as the reference normal).
                let cos = rn.dot(&face_normal);
                if (cos.abs() - 1.0).abs() > angle_tol {
                    return None;
                }
                // Project this face's origin onto the reference
                // plane's signed-distance axis. If the offset is
                // larger than `offset_tol`, the face sits on a
                // different plane and we bail out.
                let delta = *origin - ref_origin_pt;
                let signed = delta.dot(&rn);
                if signed.abs() > offset_tol {
                    return None;
                }
            }
            _ => unreachable!(),
        }
    }

    Some((ref_origin?, ref_normal?))
}

/// Detect whether `sheet`'s vertex positions lie within `max_deviation`
/// of a single best-fit plane, even when the face surface types are
/// non-planar (NURBS, sphere cap, cylinder section).
///
/// Returns `Some((centroid, normal))` when the largest distance from
/// any vertex to the best-fit plane is at most `max_deviation`. Returns
/// `None` when the sheet has fewer than 3 non-collinear vertices, or
/// when any vertex exceeds the deviation envelope.
///
/// This is the relaxed fallback for [`is_planar_sheet`]: the strict
/// path only matches exact `SurfaceGeometry::Plane` faces, while this
/// relaxed path catches sheets whose vertices happen to be co-planar
/// regardless of the underlying surface type. Used by the curved-sheet
/// multi-piece dispatcher to upgrade "effectively planar" sheets to
/// the half-space split path.
///
/// Implementation: least-squares plane fit via the covariance matrix's
/// smallest eigenvector. We compute the 3×3 covariance of the vertex
/// positions relative to the centroid, then approximate the smallest
/// eigenvector by the cross product of the two longest spans — good
/// enough for "are all these points approximately co-planar?" without
/// needing a full SVD.
pub fn is_approximately_planar_sheet(
    sheet: &BRep,
    max_deviation: f64,
) -> Option<(popkern_core::geometry::Point3, popkern_core::geometry::Vec3)> {
    use popkern_core::geometry::{Point3, Vec3};

    // Collect vertex positions. Need at least 3 to define a plane.
    let positions: Vec<Point3> = sheet.vertices.values().map(|v| v.position).collect();
    if positions.len() < 3 {
        return None;
    }

    // Centroid.
    let n = positions.len() as f64;
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for p in &positions {
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }
    let centroid = Point3::new(cx / n, cy / n, cz / n);

    // Find the vertex farthest from the centroid — defines the first
    // in-plane axis.
    let mut axis_a: Option<Vec3> = None;
    let mut best_dist_sq = 0.0;
    for p in &positions {
        let d = *p - centroid;
        let ds = d.length_squared();
        if ds > best_dist_sq {
            best_dist_sq = ds;
            axis_a = Some(d);
        }
    }
    let axis_a = axis_a?.normalized()?;

    // Find the vertex that's farthest from the line (centroid, axis_a)
    // — this gives the second axis perpendicular to axis_a.
    let mut axis_b: Option<Vec3> = None;
    let mut best_perp_dist_sq = 0.0;
    for p in &positions {
        let d = *p - centroid;
        let projected = axis_a * d.dot(&axis_a);
        let perp = d - projected;
        let perp_ds = perp.length_squared();
        if perp_ds > best_perp_dist_sq {
            best_perp_dist_sq = perp_ds;
            axis_b = Some(perp);
        }
    }
    let axis_b = axis_b?.normalized()?;

    // Plane normal = axis_a × axis_b.
    let normal = axis_a.cross(&axis_b).normalized()?;

    // Verify every vertex lies within `max_deviation` of the plane.
    for p in &positions {
        let d = *p - centroid;
        let signed = d.dot(&normal);
        if signed.abs() > max_deviation {
            return None;
        }
    }

    // Pick an orientation for the normal: if any face's surface
    // provides an outward direction, align our normal with it so the
    // downstream half-space split uses a meaningful outward direction.
    // For sheets with no usable surface data, we keep the arbitrary
    // cross-product direction — the split still produces two halves,
    // they're just labeled in an arbitrary order.
    let oriented_normal = {
        let mut chosen = normal;
        for face in sheet.faces.values() {
            let Some(surf_on_face) = face.surface.as_ref() else {
                continue;
            };
            use popkern_core::geometry::surface::SurfaceEval;
            // Sample the surface's normal at the midpoint of its
            // valid parameter range.
            let ((u_min, u_max), (v_min, v_max)) = surf_on_face.geometry.parameter_range();
            let u_mid = if u_min.is_finite() && u_max.is_finite() {
                (u_min + u_max) * 0.5
            } else {
                0.0
            };
            let v_mid = if v_min.is_finite() && v_max.is_finite() {
                (v_min + v_max) * 0.5
            } else {
                0.0
            };
            let face_n_raw = surf_on_face.geometry.normal_at(u_mid, v_mid);
            let face_n = if surf_on_face.orientation {
                face_n_raw
            } else {
                face_n_raw * -1.0
            };
            if let Some(face_n) = face_n.normalized()
                && chosen.dot(&face_n) < 0.0
            {
                chosen = chosen * -1.0;
            }
            break;
        }
        chosen
    };

    Some((centroid, oriented_normal))
}
