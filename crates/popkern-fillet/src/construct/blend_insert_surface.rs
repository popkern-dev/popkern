//! Surface assignment and reporting helpers for vertex blend insertion.

use super::*;

/// Phase 4: Assign the blend surface to the merged face.
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
/// loop); non-uniform → Gregory cap path (substitution + Path A
/// vertex relocation).
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

/// Sphere-owned corner pre-screen (PK-222): uniform radii with
/// all-analytic strips always serve through the rolling-ball sphere,
/// so the N=3 N-form dispatch skips its attempt there outright — the
/// attempted-then-declined sector complex cost a 12x allocation
/// regression on the 12-edge cube (§18 profile).
pub(super) fn sphere_owns_corner(
    radii_uniform: bool,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
) -> bool {
    use popkern_core::geometry::surface::SurfaceGeometry;
    radii_uniform
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
