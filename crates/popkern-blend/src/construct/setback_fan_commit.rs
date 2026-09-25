//! The setback cap fan and its commit postconditions, shared by the
//! corner-setback carve and the saddle setback ring (PK-527).
//! Extracted from `setback_vertex_topology` for the file budget.

use super::*;

/// Build the KP20 cap fan over `blend_face` (one ring side per edge of
/// `side_edges`, in the border's side order), author its pcurves, and
/// run the commit postconditions: seam quality, folds, the strip
/// conflict class, and seam G2 under a G2 request. On success the fan
/// is recorded with its boundary tolerances.
#[allow(clippy::too_many_arguments)]
pub(super) fn commit_cap_fan(
    brep: &mut BRep,
    blend_face: FaceHandle,
    side_edges: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
    presolved: Option<&crate::construct::setback_cap_presolve::SetbackCapSolve>,
    gate_seam_g2: bool,
) -> Result<Vec<FaceHandle>, String> {
    // ---- cap face fan ---------------------------------------------------
    let cap_faces = {
        let faces = build_cap_face_fan(
            brep,
            blend_face,
            side_edges,
            fillet_surfaces,
            vertex_position,
            presolved,
        )?;
        // Rim refit: exact cap-row edge curves (PK-096).
        refit_cap_rim(brep, &faces)?;
        faces
    };

    // ---- pcurves for the new topology -----------------------------------
    // Plane faces: rebuild loop pcurves from the 3D curves (the
    // corner path edges live on the lateral planes). Cap pieces and
    // sweep faces: the new edges follow iso/straight uv paths, so
    // endpoint-projected Linear pcurves are exact.
    author_setback_pcurves(brep, &cap_faces)?;

    // Commit postcondition: the fan's live seams must be commercial.
    // The all-or-nothing contract belongs HERE, on measured seam
    // quality, not only on structural validity — a blunt apex (18°
    // half-angle, pentagon5_blunt) builds a structurally valid fan
    // whose wide flat sectors crease 2.9e-1 at one seam corner, and
    // declining lets the sphere-fit rung serve the gentle
    // corner instead. Sharp-apex fixtures measure g0 ≤ 2e-10 and
    // g1 ≤ 2e-3 here — orders under the gates.
    let (seam_g0, seam_g1, seam_g2) =
        setback_fan_seam_refresh::worst_fan_seam_quality(brep, &cap_faces);
    if popkern_core::trace::trace_enabled("setback") {
        eprintln!("  [seam-postcondition] g0={seam_g0:.3e} g1={seam_g1:.3e} g2={seam_g2:.3e}");
    }
    if seam_g0 > 1.0e-5 || seam_g1 > 5.0e-2 {
        return Err(format!(
            "fan seam postcondition failed: g0={seam_g0:.3e} (gate 1e-5) \
             g1={seam_g1:.3e} (gate 5e-2)"
        ));
    }
    // The cap's boundary with the strips and supports is designed G1
    // too, and the cap-cap measure above never sees it. A saddle bridge
    // regime whose cap-cap seams were exact creased 7° at the cap
    // boundary and served (PK-498): gate the rim's normal angle as well.
    let rim_g1 = worst_cap_boundary_g1(brep, &cap_faces);
    if popkern_core::trace::trace_enabled("setback") {
        eprintln!("  [rim-postcondition] g1={rim_g1:.3e}");
    }
    if rim_g1 > 5.0e-2 {
        return Err(format!(
            "fan rim postcondition failed: cap-boundary g1={rim_g1:.3e} (gate 5e-2)"
        ));
    }
    // Fold postcondition: the fast structural validate() the commit
    // runs excludes self-intersection, so a fan whose pieces FOLD
    // would serve and only die later at the Part serving gate
    // (measured on the N=3 cube corner: 3 near-pole pieces fold with
    // interior normal reversals to dot -1.0 while the assembly's
    // unsigned min-regularity 8.7e-5 slid past a gate inflated by the
    // creased-border corner mismatch). Decline HERE with the class
    // named; the vertex-blend path serves the corner instead.
    // Read-only per-piece measurement; ordered collect keeps the
    // first-folded-face report identical to the serial walk (PK-302).
    let fold_of = |fh: &FaceHandle| -> bool {
        brep.faces
            .get(*fh)
            .and_then(|f| f.surface.as_ref())
            .is_some_and(|sof| match &sof.geometry {
                popkern_core::geometry::SurfaceGeometry::Gregory(patch) => {
                    patch.self_intersects_with_tolerance(brep.tolerance.linear().max(1e-6))
                }
                _ => false,
            })
    };
    #[cfg(feature = "parallel")]
    let folds: Vec<bool> = {
        use rayon::prelude::*;
        cap_faces.par_iter().map(fold_of).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let folds: Vec<bool> = cap_faces.iter().map(fold_of).collect();
    for (&fh, folded) in cap_faces.iter().zip(folds) {
        if folded {
            return Err(format!(
                "fan fold postcondition failed: cap piece {fh:?} self-intersects \
                 (assembly-born fold — creased-border wide-sector class)"
            ));
        }
    }
    // Strip-conflict postcondition (PK-225): the setback extends the
    // strips to the arcs, and flat conic sections at small radii
    // genuinely interpenetrate near the shared contact corner (see
    // setback_commit_gate module docs). The flat-section class is a
    // TRIGGER, not the verdict (PK-518): every G2 strip is a little
    // flatter than circular by construction (measured rho 0.4351 vs
    // 0.4443 at the r = 0.16 triangular-pyramid apex) and serves
    // self-intersection-clean, while the conic rho = 0.35 cube
    // corner interpenetrates. When the class fires, the committed
    // corner region is measured for self-intersection and only a
    // real intersection declines.
    if let Some(conflict) = crate::construct::setback_commit_gate::flat_section_class(fillet_surfaces)
        && let Err(intersection) = crate::apex_construction::exact_boundary_cap::
            cap_quality_postcondition::validate_candidate_region_self_intersection(
                brep,
                &cap_faces,
            )
    {
        return Err(format!(
            "strip-conflict postcondition failed: {conflict}; measured: {intersection}"
        ));
    }
    // A G2 REQUEST additionally gates the seam's normal-curvature
    // residual (PK-122): a blunt apex (18° half-angle, pentagon5_blunt)
    // builds a G0/G1-commercial fan whose wide flat sectors still jump
    // 7.0e-2 in curvature across seams; committing it would report an
    // achieved G2 the geometry does not have. Declining hands the
    // gentle corner to the sphere-fit rung. Sharp-apex fixtures
    // measure well inside the gate (the kp_polar commercial gates
    // assert the same 5e-2 bound on their live seams).
    if gate_seam_g2 && seam_g2 > 5.0e-2 {
        return Err(format!(
            "fan seam G2 postcondition failed under a G2 request: \
             g2={seam_g2:.3e} (gate 5e-2, g0={seam_g0:.3e} g1={seam_g1:.3e})"
        ));
    }

    // Record the fit residual the rim reconciliation legitimately
    // left as per-edge tolerances so the serving curve-surface
    // validation (check 9) can accept the fan — see
    // `record_boundary_tolerances` for the contract and for why this
    // must never feed the tessellator.
    setback_boundary_tolerance::record_boundary_tolerances(brep, &cap_faces);

    Ok(cap_faces)
}

/// Samples per rim edge for [`worst_cap_boundary_g1`].
const RIM_SAMPLES: usize = 33;

/// Worst sine of the normal angle across the cap's boundary: every edge
/// between a cap piece and a non-cap face, sampled pcurve-true on both
/// sides. An edge that cannot be measured counts as a full crease, so
/// it can never pass silently.
fn worst_cap_boundary_g1(brep: &BRep, cap_faces: &[FaceHandle]) -> f64 {
    use crate::cap_seam_measure::{edge_seam_faces, measure_edge_seam, resolve_edge_seam};
    let caps: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let mut worst = 0.0_f64;
    for (_, edge) in brep.edges.iter() {
        let Some(adjacency) = edge_seam_faces(brep, edge) else {
            continue;
        };
        if caps.contains(&adjacency.face_a) == caps.contains(&adjacency.face_b) {
            continue;
        }
        worst = worst.max(match resolve_edge_seam(brep, edge, &adjacency) {
            Some(seam) => measure_edge_seam(&seam, RIM_SAMPLES).g1,
            None => 1.0,
        });
    }
    worst
}
