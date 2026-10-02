//! Surface assignment and reporting helpers for vertex blend insertion.

use super::*;

/// Assign the blend surface to the merged face.
///
/// R8 C6 — takes a pre-computed `blend_boundary` because V may
/// already have been `kev`'d by `merge_corner_triangles`.
/// `build_blend_surface` only reads `boundary.vertex_position`, so
/// the removed handle in `boundary.vertex` is diagnostic only.
///
/// Returns `(method, achieved_continuity)` from the construction.
pub(super) fn assign_blend_surface(
    brep: &mut BRep,
    blend_face: FaceHandle,
    blend_boundary: &blend::BlendBoundary,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    classification: crate::vertex_classifier::VertexBlendClassification,
    tolerance: &Tolerance,
) -> Result<(CornerMethod, ContinuityLevel), String> {
    let blend_result = blend::build_blend_surface(
        blend_boundary,
        fillet_surfaces,
        brep,
        VertexMode::Smooth,
        classification,
        tolerance,
    )
    .map_err(|e| format!("blend surface computation failed: {e}"))?;

    let method = blend_result.method;
    let continuity = blend_result.continuity;
    brep.faces[blend_face].surface = Some(SurfaceOnFace {
        geometry: blend_result.surface,
        orientation: true,
    });

    Ok((method, continuity))
}

/// Minimum incident-fillet count for the setback rung at this corner
/// (`usize::MAX` disables it). Single source of truth shared by the
/// serial dispatch in `insert_vertex_blend` and the PK-278 parallel
/// presolve — the two must agree or the presolve wastes work /
/// misses corners.
pub(crate) fn setback_min_n_for(
    edge_indices: &[usize],
    edge_geometries: &[crate::construct::EdgeGeometry],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    classification: crate::vertex_classifier::VertexBlendClassification,
    tolerance: &Tolerance,
    explicit_subpatch_family: bool,
) -> usize {
    if crate::dispatch_overrides::SETBACK_DISABLE.get() == Some(true) {
        usize::MAX
    } else if !sphere_owns_corner(
        corner_radii_uniform(edge_indices, edge_geometries, tolerance),
        fillet_surfaces,
        classification,
    ) || (!explicit_subpatch_family
        && oblique_cylinder_corner(edge_indices, edge_geometries, fillet_surfaces))
    {
        3
    } else {
        5
    }
}

/// Whether a uniform cylindrical corner's strip axes are NOT mutually
/// orthogonal (PK-597). The rolling-ball sphere is only constructed for
/// the orthogonal (cube-class) trihedron: at an oblique corner — every
/// regular pyramid apex — `try_spherical_blend`'s fit on co-circular
/// samples is degenerate (PK-229/233), so a "sphere-owned" oblique
/// corner fell to the boundary-ring Gregory cap (n = 4 G1: 27° crease,
/// body G0). Such corners take the setback rung instead (measured: n = 4
/// h = 1.2 r = 0.16 G1 certifies G1 at g1 9.3e-10). A sub-patch family
/// flag keeps these corners on the vertex-blend ladder, whose reaching
/// fixture they are (direct-apex promotion stands down likewise).
fn oblique_cylinder_corner(
    edge_indices: &[usize],
    edge_geometries: &[crate::construct::EdgeGeometry],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
) -> bool {
    use popkern_core::geometry::surface::SurfaceGeometry;
    fillet_surfaces
        .iter()
        .all(|fs| matches!(fs.surface, SurfaceGeometry::Cylinder { .. }))
        && !crate::construct::direct_apex_promotion::cylinder_axes_mutually_orthogonal(
            edge_indices,
            edge_geometries,
        )
}

/// §6c-X.A.7 helper: detect whether all incident fillets at this
/// apex have the same effective radius within tolerance.
///
/// "Uniform" here means **every** edge's `radius_at_start` and
/// `radius_at_end` are equal to a single value within
/// `tolerance.linear()`.  Variable-radius fillets (start ≠ end)
/// or any per-edge mismatch makes the corner non-uniform.
///
/// Used by the cap-kind dispatch: uniform → Sphere cap path
/// (suppress sphere-fit substitution, run topology rebuild post-
/// loop); non-uniform → Gregory cap path (substitution + apex
/// topology-rebuild vertex relocation).
pub(super) fn corner_radii_uniform(
    edge_indices: &[usize],
    edge_geometries: &[crate::construct::EdgeGeometry],
    tolerance: &Tolerance,
) -> bool {
    let tol = tolerance.linear();
    let mut anchor: Option<f64> = None;
    for &i in edge_indices {
        let Some(eg) = edge_geometries.get(i) else {
            return false;
        };
        if (eg.radius_at_start - eg.radius_at_end).abs() > tol {
            return false;
        }
        let r = eg.radius_at_start;
        match anchor {
            None => anchor = Some(r),
            Some(a) if (a - r).abs() > tol => return false,
            _ => {}
        }
    }
    anchor.is_some()
}

/// The corner's smallest fillet radius over both ends of every strip,
/// the fairness measure's unit (PK-612).
pub(super) fn corner_min_radius(
    edge_indices: &[usize],
    edge_geometries: &[crate::construct::EdgeGeometry],
) -> f64 {
    edge_indices
        .iter()
        .filter_map(|&i| edge_geometries.get(i))
        .map(|eg| eg.radius_at_start.min(eg.radius_at_end))
        .fold(f64::INFINITY, f64::min)
}

/// Sphere-owned corner pre-screen (PK-222): uniform radii with
/// all-analytic strips always serve through the rolling-ball sphere,
/// so the N=3 N-form dispatch skips its attempt there outright — the
/// attempted-then-declined sector complex cost a 12x allocation
/// regression on the 12-edge cube (§18 profile).
///
/// Only a SAME-SIDE corner (all convex or all concave) is a ball's: at
/// a saddle the ball cannot sit on one side of every support. The
/// screen once ignored convexity, so an equal-radius saddle (the step
/// with two saddles, PK-583) never reached the saddle setback at a G1
/// request and served 0.43 creases from the rungs below it.
pub(super) fn sphere_owns_corner(
    radii_uniform: bool,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    classification: crate::vertex_classifier::VertexBlendClassification,
) -> bool {
    use crate::vertex_classifier::VertexBlendClassification as VC;
    use popkern_core::geometry::surface::SurfaceGeometry;
    radii_uniform
        && matches!(classification, VC::SameSideConvex | VC::SameSideConcave)
        && fillet_surfaces.iter().all(|fs| {
            matches!(
                fs.surface,
                SurfaceGeometry::Cylinder { .. } | SurfaceGeometry::Torus { .. }
            )
        })
}

/// Reclassify `CornerPatchInfo::method` from each patch's final
/// face surface type, because `surface_fix::correct_face_surfaces`
/// can override `build_blend_surface`'s initial choice (e.g.
/// promoting a Mitred fallback plane to a fitted Sphere).
pub(crate) fn reclassify_corner_methods(
    brep: &BRep,
    corner_patches: &mut [crate::types::CornerPatchInfo],
) {
    use popkern_core::geometry::surface::SurfaceGeometry;
    for info in corner_patches.iter_mut() {
        // PnS bi-4 caps produce NURBS surfaces too, but they're
        // a distinct construction with verified G¹ — don't let the
        // generic Nurbs → Gregory inference downgrade the label.
        if info.method == CornerMethod::PolyhedralNetSpline {
            continue;
        }
        // Exact-boundary caps commit NURBS sector faces with their
        // own gated continuity claim (PK-119) — the generic
        // Nurbs → Gregory inference would erase the Coons provenance.
        if info.exact_boundary_cap {
            continue;
        }
        // Variational-fit sectors (PK-574) are NURBS faces whose level
        // comes from certificates; the inference would relabel them as
        // an uncertified G1 Gregory cap.
        if info.method == CornerMethod::PlateMin {
            continue;
        }
        let Some(sof) = brep.faces.get(info.face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        // Infer the method from the final face surface. Only
        // overwrite the method/continuity fields if the surface
        // type mismatches the original construction — that means
        // `surface_fix` promoted/demoted the surface post-hoc
        // (e.g. Mitred plane → RollingBallSphere). If the surface
        // type matches, keep the construction-reported values —
        // the §5.3 C4 Gregory path reports G1 or G2 per-case based
        // on radius uniformity, and we must not flatten that to G1.
        let inferred = match &sof.geometry {
            SurfaceGeometry::Sphere { .. } => {
                Some((CornerMethod::RollingBallSphere, ContinuityLevel::G1))
            }
            SurfaceGeometry::Nurbs(_) => Some((CornerMethod::Gregory, ContinuityLevel::G1)),
            SurfaceGeometry::Plane { .. } => Some((CornerMethod::Mitred, ContinuityLevel::G0)),
            _ => None,
        };
        if let Some((inferred_method, inferred_cont)) = inferred
            && info.method != inferred_method
        {
            info.method = inferred_method;
            info.achieved_continuity = inferred_cont;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_vertex_blend_too_few_anchors() {
        // A cube vertex has 3 neighbors, but if we pass no edge indices
        // the anchor set will be non-empty (all neighbors). To trigger the
        // "fewer than 2 anchors" path we need a degenerate vertex.
        // Instead, test with a default (non-existent) vertex handle.
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let bogus = popkern_core::topology::VertexHandle::default();
        let result = insert_vertex_blend(
            &mut brep,
            bogus,
            &[],
            &[],
            &Tolerance::default(),
            crate::types::ContinuityLevel::G1,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            None,
            None,
            &crate::construct::blend_insert_setback_user::CornerRequest::cap_only(
                crate::types::CornerCapPreference::default(),
            ),
            false,
        );
        assert!(result.is_err(), "should fail for invalid vertex");
    }

    #[test]
    fn insert_vertex_blend_on_cube_vertex() {
        // Calling insert_vertex_blend on a real cube vertex with no
        // edge geometries should fail during blend surface computation
        // (since there are no fillet surfaces to blend), but it should
        // get past the anchor check (3 anchors >= 2).
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let vertex = brep.vertices.keys().next().unwrap();
        let result = insert_vertex_blend(
            &mut brep,
            vertex,
            &[],
            &[],
            &Tolerance::default(),
            crate::types::ContinuityLevel::G1,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            None,
            None,
            &crate::construct::blend_insert_setback_user::CornerRequest::cap_only(
                crate::types::CornerCapPreference::default(),
            ),
            false,
        );
        // This should fail (no edge geometries -> empty fillet surfaces
        // -> blend boundary computation fails), but NOT at the anchor check.
        assert!(result.is_err());
        let err = result.unwrap_err();
        // The error should be from blend computation, not "fewer than 2 anchors".
        assert!(
            !err.contains("fewer than 2"),
            "should not fail at anchor check: {err}"
        );
    }
}
