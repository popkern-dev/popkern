//! The setback cap fan and its commit postconditions, shared by the
//! corner-setback carve and the saddle setback ring (PK-527).
//! Extracted from `setback_vertex_topology` for the file budget.

use super::*;

/// Build the KP20 cap fan over `blend_face` (one ring side per edge of
/// `side_edges`, in the border's side order), author its pcurves, and
/// run the commit postconditions: seam quality, folds and the strip
/// conflict class. On success the fan is recorded with its boundary
/// tolerances; its seam-G2 verdict is [`fan_seam_g2_failure`].
#[allow(clippy::too_many_arguments)]
pub(super) fn commit_cap_fan(
    brep: &mut BRep,
    blend_face: FaceHandle,
    side_edges: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
    presolved: Option<&crate::construct::setback_cap_presolve::SetbackCapSolve>,
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
    author_setback_pcurves(brep, &cap_faces, &[])?;

    // Commit postcondition: the fan's live seams must be commercial.
    // The all-or-nothing contract belongs HERE, on measured seam
    // quality, not only on structural validity — a blunt apex (18°
    // half-angle, pentagon5_blunt) builds a structurally valid fan
    // whose wide flat sectors crease 2.9e-1 at one seam corner, and
    // declining lets the sphere-fit rung serve the gentle
    // corner instead. Sharp-apex fixtures measure g0 ≤ 2e-10 and
    // g1 ≤ 2e-3 here — orders under the gates.
    let (seam_g0, seam_g1, _) = setback_fan_seam_refresh::worst_fan_seam_quality(brep, &cap_faces);
    if popkern_core::trace::trace_enabled("setback") {
        eprintln!("  [seam-postcondition] g0={seam_g0:.3e} g1={seam_g1:.3e}");
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
    let jets = crate::cap_seam_side::SurfaceJetCache::default();
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
        worst = worst.max(match resolve_edge_seam(brep, edge, &adjacency, &jets) {
            Some(seam) => {
                measure_edge_seam(
                    &seam,
                    RIM_SAMPLES,
                    crate::cap_seam_continuity::SEAM_END_MARGIN,
                )
                .g1
            }
            None => 1.0,
        });
    }
    worst
}

/// The committed fan's seam-G2 verdict under a G2 request:
/// `Some(reason)` when it reads below G2, and the rung then serves the
/// same fan as a measured G1 downgrade (no rebuild: the fan does not
/// depend on the request).
///
/// Two instruments, either failing fails the fan:
///
/// * the fan's own dense cap-cap sampling (PK-122): a blunt apex (18°
///   half-angle, pentagon5_blunt) builds a G0/G1-commercial fan whose
///   wide flat sectors jump 7.0e-2 in curvature across seams in a
///   narrow spike;
/// * the whole blend's continuity evidence over the fan's cap-cap and
///   rim seams (PK-596): a pyramid-4 fan at r = 0.05 read under the
///   gate on the first instrument (end margin, finite-difference
///   normals) while its spokes certifiably jump 8.7e-2 near the
///   centre, and the corner reported a G2 the body did not have.
///
/// The evidence is sampled first, then certified where it matters:
/// every gated seam whose sampled reading fails is certified (PK-499,
/// decision mode), and its certificate confirms (`lower > gate`) or
/// overturns (`upper <= gate`) the failure exactly as the whole-blend
/// evidence decides it. Seams sampling under the gate are not
/// certified here (the fan's 40 seams certify in ~2.3 s against ~1 s
/// for the whole corner); the served body's evidence certifies them.
pub(crate) fn fan_seam_g2_failure(brep: &BRep, cap_faces: &[FaceHandle]) -> Option<String> {
    let (_, _, seam_g2) = setback_fan_seam_refresh::worst_fan_seam_quality(brep, cap_faces);
    if seam_g2 > 5.0e-2 {
        return Some(format!("fan seam sampling g2={seam_g2:.3e} (gate 5e-2)"));
    }

    use crate::ContinuityLevel;
    use crate::blend_continuity_certification::{
        ContinuityCertification, certify_seam_for_evidence, evidence_certify_options,
    };
    use crate::blend_continuity_evidence::{
        ContinuityGates, measure_blend_continuity_evidence_with,
    };
    let gates = ContinuityGates::commercial();
    let mut evidence = measure_blend_continuity_evidence_with(
        brep,
        cap_faces,
        ContinuityLevel::G2,
        &gates,
        RIM_SAMPLES,
        ContinuityCertification::Deferred,
    );
    if evidence.achieved >= ContinuityLevel::G2 {
        return None;
    }
    let options = evidence_certify_options(brep, &gates);
    for seam in evidence.seams.iter_mut().filter(|s| s.gated) {
        let (g0, g1, g2) = seam.classified_residuals(&gates);
        if g0 > gates.g0_gap || g1 > gates.g1_normal || g2 > gates.g2_curvature {
            seam.certificate = Some(certify_seam_for_evidence(brep, seam.edge, &options));
            seam.certification_deferred = false;
        }
    }
    evidence.summarize();
    if evidence.achieved >= ContinuityLevel::G2 {
        return None;
    }
    let (seam, (g0, g1, g2)) = evidence
        .seams
        .iter()
        .filter(|s| s.gated)
        .map(|s| (s, s.classified_residuals(&gates)))
        .max_by(|a, b| a.1.2.total_cmp(&b.1.2))?;
    Some(format!(
        "{} {:?} (seam {:?} {:?}: g0={g0:.3e} g1={g1:.3e} g2={g2:.3e}, gate 5e-2)",
        if seam.certificate_decides(&gates) {
            "certified"
        } else {
            "sampled"
        },
        evidence.achieved,
        seam.edge,
        seam.kind
    ))
}
