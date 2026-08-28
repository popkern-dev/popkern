//! Transactional PnS dispatch acceptance for vertex-blend insertion.

use popkern_core::error::ValidationError;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend::BlendBoundary;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo};
use crate::vertex_classifier::VertexBlendClassification;

use super::blend_insert_diagnostics::pns_dispatch_declined_diagnostic;

pub(super) struct PnsDispatchAcceptanceInput<'a> {
    pub brep: &'a mut BRep,
    pub blend_face: FaceHandle,
    pub boundary: &'a BlendBoundary,
    pub fillet_surfaces: &'a [FilletSurfaceResult],
    pub vertex: VertexHandle,
    pub classification: VertexBlendClassification,
    pub n_sides: usize,
    pub diagnostics: &'a mut Vec<ApexConstructionDiagnostic>,
}

pub(super) fn try_pns_dispatch_with_acceptance(
    input: PnsDispatchAcceptanceInput<'_>,
) -> Option<Vec<CornerPatchInfo>> {
    let PnsDispatchAcceptanceInput {
        brep,
        blend_face,
        boundary,
        fillet_surfaces,
        vertex,
        classification,
        n_sides,
        diagnostics,
    } = input;
    let pns_snapshot = brep.clone();
    let patches = match crate::blend::pns_promote::try_pns_dispatch(
        brep,
        blend_face,
        boundary,
        fillet_surfaces,
        vertex,
        classification,
        boundary.bridge_fallbacks.clone(),
    ) {
        Ok(patches) => patches,
        Err(error) => {
            // Restore before returning: a partially-inserted cap must
            // not survive a decline. Then record WHY — this path used
            // to return None silently, so every construction failure
            // looked identical to "PnS was never tried".
            *brep = pns_snapshot;
            diagnostics.push(pns_dispatch_declined_diagnostic(
                n_sides,
                "construction_failed",
                Some(error.to_string()),
            ));
            return None;
        }
    };

    let cap_faces = patches.iter().map(|patch| patch.face).collect::<Vec<_>>();
    // Measure pcurve fidelity on the cap faces and their neighbours —
    // the faces the acceptance gate judges — with the validated
    // instrument. Read-only; reported through the non-blocking
    // channel so it cannot change which route serves the corner.
    let mut audited: Vec<popkern_core::topology::FaceHandle> = cap_faces.clone();
    for &cap in &cap_faces {
        let Some(face_data) = brep.faces.get(cap) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut he = start;
        for _ in 0..32 {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(twin) = brep.half_edges.get(hd.twin)
                && let Some(l) = brep.loops.get(twin.loop_)
                && !audited.contains(&l.face)
            {
                audited.push(l.face);
            }
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    // PK-192: the acceptance gate below judges the cap against its
    // neighbours, but PnS dispatch runs BEFORE the pipeline's
    // pcurve-completion passes — so those neighbours are unfinished at
    // this moment. Measured on the same class of face: 3.277 here
    // against 4.140e-7 at pipeline end, a factor of eight million. The
    // cap was being refused for its neighbours' unfinished state.
    //
    // Complete them first, so the gate measures geometry rather than
    // timing.
    // Verified: this succeeds on every audited face, including the one
    // the acceptance gate then reports against — so a silent recompute
    // failure is NOT the explanation for the two measurements
    // disagreeing (see the open inconsistency on PK-192).
    let mut recompute_failures = 0usize;
    // Snapshot each audited face's half-edge pcurve images so we can
    // say whether the recompute actually CHANGED anything, rather than
    // succeeding while reproducing what was already there.
    let before: Vec<(
        popkern_core::topology::HalfEdgeHandle,
        popkern_core::geometry::Point3,
    )> = audited
        .iter()
        .filter_map(|&f| brep.faces.get(f))
        .filter_map(|fd| brep.loops.get(fd.outer_loop))
        .flat_map(|ld| {
            let mut out = Vec::new();
            let start = ld.half_edge;
            let mut he = start;
            for _ in 0..32 {
                let Some(hd) = brep.half_edges.get(he) else {
                    break;
                };
                if let Some(pc) = hd.pcurve.as_ref() {
                    out.push((he, pc.point_at(0.5 * (pc.t_start + pc.t_end))));
                }
                he = hd.next;
                if he == start {
                    break;
                }
            }
            out
        })
        .collect();
    for &face in &audited {
        // At the low valences the cap faces carry analytically-
        // authored, rate-matched pcurves (see
        // `author_pns_cap_edge_geometry`); the generic 2-D projection
        // recompute lands interior preimages on bi-4 patches and must
        // not overwrite them. Neighbours still need completing — that
        // is the timing defect this pass exists for. The n >= 5 route
        // keeps its long-standing recompute.
        if n_sides <= 4 && cap_faces.contains(&face) {
            continue;
        }
        if brep.recompute_face_pcurves(face).is_err() {
            recompute_failures += 1;
        }
    }
    let mut changed = 0usize;
    for (he, mid) in &before {
        if let Some(pc) = brep.half_edges.get(*he).and_then(|h| h.pcurve.as_ref()) {
            let now = pc.point_at(0.5 * (pc.t_start + pc.t_end));
            if (now.x - mid.x).abs() > 1e-12 || (now.y - mid.y).abs() > 1e-12 {
                changed += 1;
            }
        }
    }
    // Direct check, no pcurve in between: does each edge's own 3D
    // curve actually LIE on the surface of the face beside it? If the
    // best pcurve the kernel can compute still sits 4.6e-2 off, the
    // curve and the surface may simply disagree — which would be a G0
    // defect underneath, not a trim-curve problem.
    //
    // Probe-gated (PK-298): this audit is PK-203 investigation
    // instrumentation, not a gate input — nothing downstream reads
    // it. Its 21 Newton `closest_point` projections per half-edge
    // across every audited face measured 42ms of the 52ms n=3
    // dispatch attempt (~80%), and with PNS_MIN_VALENCE=3 that cost
    // lands on every cube-class corner that ultimately declines —
    // the 12-edge cube fillet regressed 106ms → 515ms past its CI
    // budget. `POPKERN_PNS_CURVE_AUDIT=1` restores the measurement;
    // the diagnostic records that the audit was skipped otherwise.
    let mut worst_curve_off_surface = 0.0_f64;
    // Identity of the worst curve/surface disagreement — the scalar
    // alone cannot say WHOSE edge is off WHICH surface, which is the
    // whole PK-203 question.
    let mut worst_off_detail = if pns_curve_audit_enabled() {
        String::from("none")
    } else {
        String::from("skipped:set POPKERN_PNS_CURVE_AUDIT=1")
    };
    let audit_faces: &[popkern_core::topology::FaceHandle] = if pns_curve_audit_enabled() {
        &audited
    } else {
        &[]
    };
    for &face in audit_faces {
        let Some(fd) = brep.faces.get(face) else {
            continue;
        };
        let Some(sf) = fd.surface.as_ref() else {
            continue;
        };
        let Some(ld) = brep.loops.get(fd.outer_loop) else {
            continue;
        };
        let start = ld.half_edge;
        let mut he = start;
        for _ in 0..32 {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(curve) = brep.edges.get(hd.edge).and_then(|e| e.curve.as_ref()) {
                for k in 0..=20 {
                    let f = 0.05 + 0.90 * f64::from(k) / 20.0;
                    let p = {
                        use popkern_core::geometry::curve::CurveEval;
                        (&curve.geometry as &dyn CurveEval)
                            .point_at(curve.t_start + f * (curve.t_end - curve.t_start))
                    };
                    use popkern_core::geometry::surface::SurfaceEval;
                    if let Some((u, v)) = sf.geometry.closest_point(&p) {
                        let d = sf.geometry.point_at(u, v).distance_to(&p);
                        if d > worst_curve_off_surface {
                            worst_curve_off_surface = d;
                            let is_cap = cap_faces.contains(&face);
                            worst_off_detail = format!(
                                "edge={:?} face={face:?}(cap={is_cap},{}) curve={} f={f:.2} at=({:.4},{:.4},{:.4})",
                                hd.edge,
                                surface_kind_name(&sf.geometry),
                                curve_kind_name(&curve.geometry),
                                p.x,
                                p.y,
                                p.z,
                            );
                        }
                    }
                }
            }
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    let fidelity = crate::pcurve_fidelity::measure_pcurve_fidelity(brep, &audited);
    // A bi-4 patch has domain [0,1]^2 and EXTRAPOLATES violently
    // outside it, so a large fidelity reading on a cap face is
    // ambiguous until the domain is checked: a pcurve straying out of
    // domain is a real defect, an instrument sampling out of domain is
    // not. Report the UV extent so the two are distinguishable.
    let mut uv_min = (f64::MAX, f64::MAX);
    let mut uv_max = (f64::MIN, f64::MIN);
    for &cap in &cap_faces {
        let Some(fd) = brep.faces.get(cap) else {
            continue;
        };
        let Some(ld) = brep.loops.get(fd.outer_loop) else {
            continue;
        };
        let start = ld.half_edge;
        let mut he = start;
        for _ in 0..32 {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(pc) = hd.pcurve.as_ref() {
                for k in 0..=10 {
                    let t = pc.t_start + f64::from(k) / 10.0 * (pc.t_end - pc.t_start);
                    let uv = pc.point_at(t);
                    uv_min = (uv_min.0.min(uv.x), uv_min.1.min(uv.y));
                    uv_max = (uv_max.0.max(uv.x), uv_max.1.max(uv.y));
                }
            }
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    diagnostics.push(pns_dispatch_declined_diagnostic(
        n_sides,
        "pcurve_fidelity_audit",
        Some(format!(
            "{fidelity:?};cap_uv_extent=({:.4},{:.4})..({:.4},{:.4});recompute_failures={};pcurves_changed={}/{}",
            uv_min.0,
            uv_min.1,
            uv_max.0,
            uv_max.1,
            recompute_failures,
            changed,
            before.len()
        ) + &format!(
            ";worst_curve_off_surface={worst_curve_off_surface:.3e}[{worst_off_detail}]"
        )),
    ));
    crate::construct::corner_patch_support_completion::complete_external_cap_support_pcurves(
        brep,
        &patches,
        diagnostics,
    );
    // The completion pass "repairs" cap-boundary pcurves with 2-D
    // projection, which is measurably worse than the analytic
    // authoring on bi-4 patches (3.1e-2 vs 6.5e-3). Its support-side
    // work is wanted; its cap-side work is re-authored. Low-valence
    // only — the n >= 5 route never had analytic cap pcurves.
    if n_sides <= 4 {
        crate::blend::pns_promote::reauthor_pns_cap_edge_geometry(brep, &cap_faces);
    }
    // PK-192: rim/spoke split of the CAP FACES ONLY, measured at the
    // same pipeline moment the acceptance gate below judges them —
    // after pcurve recompute and support completion, before rollback.
    // The all-half-edge audit above spans cap AND neighbour faces, so
    // its worst reading cannot attribute a deviation to the cap's own
    // border (rim) versus its internal seams (spokes); this one can.
    // Measured on the final body instead, these numbers describe the
    // incumbent cap, not PnS — that mistake is why this exists.
    let rim_spoke = crate::pcurve_fidelity::measure_cap_boundary_fidelity_split(brep, &cap_faces);
    let describe = |family: Option<crate::pcurve_fidelity::PcurveFidelity>| {
        family
            .map(|f| {
                format!(
                    "worst={:.4e} nearest={:.4e} edge={:?} face={:?} fraction={:.2} measured={}",
                    f.worst,
                    f.worst_nearest_point,
                    f.worst_edge,
                    f.worst_face,
                    f.worst_fraction,
                    f.measured
                )
            })
            .unwrap_or_else(|| "unmeasured".to_string())
    };
    diagnostics.push(pns_dispatch_declined_diagnostic(
        n_sides,
        "cap_rim_spoke_split",
        Some(format!(
            "rim[{}];spoke[{}]",
            describe(rim_spoke.rim),
            describe(rim_spoke.spoke)
        )),
    ));
    // Structural variant (PK-273): the filter below only matches
    // structural classes — the full variant's Gregory sweep is waste.
    let topology_error = brep
        .validate_collect_structural()
        .into_iter()
        .find(|error| {
            matches!(
                error,
                ValidationError::BrokenCycle { .. }
                    | ValidationError::ShellGap { .. }
                    | ValidationError::OpenLoop { .. }
            )
        });
    let continuity_error = crate::construct::cap_outer_acceptance::cap_support_acceptance_error(
        brep, &cap_faces, 5.0e-3, 5.0e-2,
    );
    if topology_error.is_none() && continuity_error.is_none() {
        return Some(patches);
    }

    if let Some(diagnostic) =
        crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic(
            brep, &cap_faces, 5.0e-3, 5.0e-2,
        )
    {
        diagnostics.push(diagnostic);
    }
    // PK-203 timing probe: every quantity this gate refuses on has
    // been measured to improve by orders of magnitude between this
    // moment and pipeline end (support 110,000x; the cap rim 3x from
    // one completion pass alone) — but the refused cap never reaches
    // those passes, so nobody has ever measured a FINISHED PnS cap.
    let detail_code = topology_error
        .map(|error| format!("topology_validation_failed={error}"))
        .or(continuity_error);
    diagnostics.push(pns_dispatch_declined_diagnostic(
        n_sides,
        "cap_support_acceptance",
        detail_code,
    ));
    *brep = pns_snapshot;
    None
}

/// Short curve-kind label for the audit's worst-sample identity.
fn curve_kind_name(curve: &popkern_core::geometry::CurveGeometry) -> &'static str {
    use popkern_core::geometry::CurveGeometry;
    match curve {
        CurveGeometry::Line { .. } => "Line",
        CurveGeometry::Arc { .. } => "Arc",
        CurveGeometry::Circle { .. } => "Circle",
        CurveGeometry::Ellipse { .. } => "Ellipse",
        CurveGeometry::Parabola { .. } => "Parabola",
        CurveGeometry::Hyperbola { .. } => "Hyperbola",
        CurveGeometry::Nurbs(_) => "Nurbs",
    }
}

/// Short surface-kind label for the audit's worst-sample identity.
fn surface_kind_name(surface: &popkern_core::geometry::SurfaceGeometry) -> &'static str {
    use popkern_core::geometry::SurfaceGeometry;
    match surface {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
        SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
    }
}

/// PK-203 diagnostic-probe switch: `POPKERN_PNS_CURVE_AUDIT=1`
/// enables the curve-off-surface audit inside the PnS acceptance
/// harness (21 Newton projections per half-edge over the audited
/// faces — ~80% of a declining n=3 attempt's cost, see PK-298).
/// Default off; the acceptance decision never reads its output.
fn pns_curve_audit_enabled() -> bool {
    std::env::var("POPKERN_PNS_CURVE_AUDIT")
        .map(|v| v == "1")
        .unwrap_or(false)
}
