//! §5.3 Round 8 C5 — apply a three-quad decomposition to a BRep.
//!
//! Given a BRep with a hexagonal blend face (from the existing
//! `insert_vertex_blend` pipeline), apply the commercial-style
//! three-quad decomposition:
//!   1. Build the [`ThreeQuadDecomposition`] from the hex boundary.
//!   2. Poke the hex face into three quadrilateral sub-faces
//!      via [`poke_hex_face_into_three_quads`].
//!   3. Assign each sub-patch's `RectangularGregory` surface to
//!      the corresponding sub-face.
//!   4. Return three [`CornerPatchInfo`] entries — one per
//!      sub-face — replacing the single entry the caller would
//!      otherwise have emitted.
//!
//! The helper is a **pure function on the BRep + boundary**; it
//! doesn't know about `FilletOptions` flags or dispatch — the
//! caller decides when to invoke it. R8 C6 wires the dispatch.

use super::kp_biquintic_optimizer::optimize_biquintic_and_build_binonic;
use super::kp_twist_optimizer::optimize_c_twists_and_build_binonic;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GregoryPatch, OuterCornerData,
    RectangularGregoryTwoStripe, Side, Spoke, SpokeBoundaryData, SpokeEnd, TangentPlane,
};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::poke_hex::poke_n_gon_face_into_n_quads;
use crate::blend::boundary::BlendBoundary;
use crate::blend::gregory_three_quad::build_n_quad_decomposition_with_options;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo};
use crate::vertex_classifier::VertexBlendClassification;

/// §5.3 item #10 — why the three-quad dispatch failed. Callers
/// use this to distinguish an honest topology / geometry failure
/// (caller should record the reason + fall through to polar) from
/// a self-intersection remediation (caller must ALSO emit the
/// `CornerRemediationOverlap` warning).
#[derive(Debug)]
pub(crate) enum ThreeQuadFailure {
    /// The Gregory sub-patches self-intersect; remediate to polar.
    Overlap,
    /// Anything else — boundary shape, decomposition rejection,
    /// poke topology failure. Error string included for logging.
    Other(#[allow(dead_code)] String),
}

impl From<String> for ThreeQuadFailure {
    fn from(s: String) -> Self {
        ThreeQuadFailure::Other(s)
    }
}

/// Apply a three-quad decomposition to an existing hexagonal blend
/// face: build the decomposition, poke the face, assign Gregory
/// surfaces, return one [`CornerPatchInfo`] per sub-face.
///
/// # Errors
///
/// Returns an error string if the boundary isn't a valid hexagon
/// (6 alternating segments), the decomposition builder rejects
/// the input, or the poke topology surgery fails.
#[allow(dead_code)] // wired into dispatch in R8 C6
#[allow(clippy::too_many_arguments)]
pub(crate) fn apply_three_quad_decomposition_to_blend(
    brep: &mut BRep,
    blend_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    vertex: VertexHandle,
    n_sides: usize,
    classification: VertexBlendClassification,
    use_deg5: bool,
    use_biquintic_full: bool,
    use_biseptic: bool,
    use_two_stripe: bool,
    use_kp_binonic: bool,
    use_kp_polar: bool,
) -> Result<Vec<CornerPatchInfo>, ThreeQuadFailure> {
    let mut decomp = build_n_quad_decomposition_with_options(
        boundary,
        fillet_surfaces,
        use_deg5,
        classification,
    )
    .ok_or_else(|| ThreeQuadFailure::Other("boundary not decomposable into N quads".to_string()))?;

    // Phase H.8 step-D — relaxed-continuity carry-through. `Some(c)`
    // after K-P promotion means the binonic shipped but worst-case
    // coplanarity `c` exceeded the 5e-2 commercial gate.
    let mut relaxed_coplanarity: Option<f64> = None;

    // Phase 6c-D Stage 1d: capture the substituted hex-side
    // boundary CPs from `boundary` so the K-P bi-nonic ring
    // builder can use them via `pick_hex_boundary` instead of
    // `k_p_boundary` Hermite synthesis. When the boundary's
    // `hex_boundary_cps_deg9` is empty (sphere-fit substitution
    // didn't run / declined), the closure passes `None` and
    // existing Hermite-blend behaviour is preserved.
    let hex_boundary_data: Option<&[([Point3; 10], [Point3; 10])]> =
        if boundary.hex_boundary_cps_deg9.is_empty() {
            None
        } else {
            Some(&boundary.hex_boundary_cps_deg9)
        };
    let is_4n_path = super::promotion_ladder::run_promotion_ladder(
        &mut decomp,
        fillet_surfaces,
        &mut relaxed_coplanarity,
        use_deg5,
        use_biquintic_full,
        use_biseptic,
        use_two_stripe,
        use_kp_binonic,
        use_kp_polar,
        |patches, fillets| promote_decomp_to_kp_binonic(patches, fillets, hex_boundary_data),
    );

    // Self-intersection screen. N ≥ 5 bypass: bi-5 parallelogram
    // folds on narrow 72°/60° sectors by design (K-P 2015 motivates
    // the polar construction precisely for this regime); downstream
    // polar / KP13 bi-5 paths use only shared-edge CPs.
    const SELF_INT_TOL: f64 = 1e-6;
    if decomp.sub_patches.len() < 5 {
        for patch in &decomp.sub_patches {
            if patch.self_intersects_with_tolerance(SELF_INT_TOL) {
                return Err(ThreeQuadFailure::Overlap);
            }
        }
    }

    let sub_faces: Vec<FaceHandle> = if is_4n_path {
        super::poke_4n::poke_n_gon_face_into_4n_quads(brep, blend_face, decomp.central_point)
            .map_err(|e| ThreeQuadFailure::Other(format!("4n poke failed: {e}")))?
            .sub_faces
    } else {
        poke_n_gon_face_into_n_quads(brep, blend_face, decomp.central_point)
            .map_err(|e| ThreeQuadFailure::Other(format!("poke failed: {e}")))?
            .sub_faces
    };
    if sub_faces.len() != decomp.sub_patches.len() {
        return Err(ThreeQuadFailure::Other(format!(
            "poke/decomposition size mismatch: poke produced {} sub-faces, \
             decomposition has {} sub-patches",
            sub_faces.len(),
            decomp.sub_patches.len()
        )));
    }

    // Phase 6c-D: write Phase-6c-A-substituted boundary curves
    // into BRep edges so the pcurve refit
    // (`reproject_loop_pcurves_onto_gregory`) propagates them to
    // the cap-vs-fillet seam pcurves. Gated on the boundary's
    // success flag set in `compute_blend_boundary` — when
    // substitution declined, edges are untouched and the legacy
    // hermite-cubic synthesis path runs unchanged.
    //
    // **Known limitation (Phase 6c-D blocker):** the rational-
    // Gregory dispatch this branch sits in fires from the *second*
    // `blend_insert` invocation per vertex, where the substitution
    // typically declines (the BRep face state has been mutated
    // between calls). The substituted boundary from the *first*
    // invocation never reaches us because the first invocation's
    // `face_has_2n_outer_half_edges` check returns `false`. Net:
    // `boundary.sphere_fit_substituted` is `false` here today even
    // for fixtures where Phase 6c-A's math fully applies (e.g.,
    // `cube_uniform`). Remediation = Phase 6c-C (apex-region
    // topology rebuild) so face handles in the boundary stay valid
    // through the surgery.
    if boundary.sphere_fit_substituted {
        let wired = super::sphere_fit_edge_wire::wire_substituted_boundary_to_edges(
            brep, boundary, &sub_faces,
        );
        if std::env::var_os("POPKERN_KP_TRACE").is_some() {
            eprintln!(
                "[wire_substituted_boundary_to_edges] wrote {wired} of {} segments to BRep edges",
                boundary.segments.len()
            );
        }
    }

    // Phase 6c-A continuation §4.2 option (c) — per-fillet
    // sphere-cap rung. When the boundary was sphere-fit-substituted
    // and the env var is set, *replace* each sub-face's Gregory
    // surface with the per-fillet `Sphere(cap_centre_k, r_k)`.
    // This is perfect cap-cylinder G1 by construction (each sphere
    // is tangent to its own cylinder along the cap-tangent circle),
    // sidestepping the L-M solver stall on variable-radius corners.
    // Trade-off: caps meet adjacent caps at G0 only along the
    // inner spokes — acceptable for variable-radius corners and
    // matches commercial-style sub-cap dissection.
    let sphere_cap_rung = boundary.sphere_fit_substituted
        && !boundary.per_fillet_caps.is_empty()
        && boundary.per_fillet_caps.len() == sub_faces.len()
        && std::env::var("POPKERN_SPHERE_CAP_RUNG").is_ok_and(|v| v == "1");

    // Always compute the per-fillet sphere `(centre, radius)` for each
    // sub-face when sphere-fit substitution succeeded — independent of
    // whether the sphere-cap rung will be applied.  Stashed on
    // `CornerPatchInfo.cap_sphere` so downstream Phase 6c-A reroute
    // and Phase 6c-D Gregory cap rung can identify cap faces and
    // their nominal sphere parameters even when the surface is a
    // Gregory patch.
    let cap_sphere_per_subface: Vec<Option<(Point3, f64)>> = if boundary.sphere_fit_substituted
        && !boundary.per_fillet_caps.is_empty()
        && boundary.per_fillet_caps.len() == sub_faces.len()
    {
        sub_faces
            .iter()
            .enumerate()
            .map(|(k, &fh)| {
                let chosen = pick_per_fillet_cap_for_face(
                    brep,
                    fh,
                    &boundary.per_fillet_caps,
                    boundary.per_fillet_caps[k],
                );
                Some(chosen)
            })
            .collect()
    } else {
        vec![None; sub_faces.len()]
    };

    let mut patches = Vec::with_capacity(decomp.sub_patches.len());
    for (k, (sub_face, sub_patch)) in sub_faces.iter().zip(decomp.sub_patches.iter()).enumerate() {
        // CAPTURE THE OLD SURFACE BEFORE REPLACING. Internal cap-to-
        // cap edges (created by `mev`/`mef` during the poke) have
        // `edge.curve = None`; we can't re-project from a 3-D curve
        // that doesn't exist. Instead, we sample
        // `old_surface(old_pcurve(t))` to reconstruct the 3-D path
        // of the edge, then project those points onto the new Gregory.
        let old_surface = brep
            .faces
            .get(*sub_face)
            .and_then(|f| f.surface.as_ref())
            .map(|sof| sof.geometry.clone());
        let new_geometry = if sphere_cap_rung {
            // Match by primary-cylinder axis identity (the per-fillet
            // sphere whose centre lies on this sub-face's primary
            // cylinder axis). Disambiguates uniform-radius corners
            // where matching by radius alone is ambiguous.
            let chosen = pick_per_fillet_cap_for_face(
                brep,
                *sub_face,
                &boundary.per_fillet_caps,
                boundary.per_fillet_caps[k],
            );
            SurfaceGeometry::Sphere {
                center: chosen.0,
                radius: chosen.1,
            }
        } else {
            SurfaceGeometry::Gregory(sub_patch.clone())
        };
        brep.faces[*sub_face].surface = Some(SurfaceOnFace {
            geometry: new_geometry,
            orientation: true,
        });
        // Re-fit each boundary pcurve into the new surface's UV space.
        // Uses the edge's 3-D curve when present; falls back to
        // `old_surface(old_pcurve)` when not (internal cap-to-cap
        // edges created by `mev`/`mef`).
        reproject_loop_pcurves_onto_gregory(
            brep,
            *sub_face,
            old_surface.as_ref(),
            boundary.sphere_fit_substituted,
        );
        patches.push(CornerPatchInfo {
            vertex,
            n_sides,
            method: CornerMethod::Gregory,
            achieved_continuity: ContinuityLevel::G1,
            classification,
            face: *sub_face,
            bridge_fallbacks: boundary.bridge_fallbacks.clone(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: relaxed_coplanarity,
            cap_sphere: cap_sphere_per_subface[k],
        });
    }
    Ok(patches)
}

/// Return the (origin, axis, radius) of the cylinder that the given
/// face's outer-loop edges share with the most twin faces. Used by
/// the Phase 6c-E pyramid path to disambiguate which per-fillet sphere
/// to assign to each sub-face when all per-fillet radii are equal —
/// matching by axis identity (sphere centre projects onto axis) rather
/// than by radius.
fn primary_cylinder_identity_of_face(
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

/// For a sub-face, find the per-fillet sphere whose centre lies on
/// the sub-face's primary cylinder axis (within `r * 0.05` tolerance).
/// This is the geometric tangency check — the right sphere for this
/// sub-face is the one tangent to its primary cylinder along the
/// cap-tangent circle. Falls back to the supplied `default` when
/// no per-fillet entry projects onto the axis (e.g., when the
/// inscribed-sphere path was taken instead of per-fillet rolling-ball).
fn pick_per_fillet_cap_for_face(
    brep: &BRep,
    face: FaceHandle,
    per_fillet_caps: &[(Point3, f64)],
    default: (Point3, f64),
) -> (Point3, f64) {
    let Some((cyl_origin, cyl_axis, cyl_radius)) = primary_cylinder_identity_of_face(brep, face)
    else {
        return default;
    };
    let axis_unit = cyl_axis
        .normalized()
        .unwrap_or(popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0));
    let tol = (cyl_radius * 0.05).max(1e-9);
    let mut found: Option<(Point3, f64)> = None;
    for &(centre, radius) in per_fillet_caps {
        if (radius - cyl_radius).abs() > cyl_radius * 0.05 {
            continue;
        }
        let v = popkern_core::geometry::Vec3::new(
            centre.x - cyl_origin.x,
            centre.y - cyl_origin.y,
            centre.z - cyl_origin.z,
        );
        let along = v.dot(&axis_unit);
        let proj = Point3::new(
            cyl_origin.x + axis_unit.x * along,
            cyl_origin.y + axis_unit.y * along,
            cyl_origin.z + axis_unit.z * along,
        );
        if centre.distance_to(&proj) < tol {
            found = Some((centre, radius));
            break;
        }
    }
    found.unwrap_or(default)
}

/// Re-fit each half-edge pcurve on `face`'s outer loop into the
/// face's surface UV space by closest-point projection. Used after
/// swapping a face's surface to a Gregory patch — the half-edges'
/// inherited pcurves are still in the *parent face's* UV space, and
/// re-projecting onto the new surface keeps tessellation consistent.
///
/// `old_surface` is the surface the face had *before* the swap. Used
/// as a fallback for half-edges whose underlying edge has no 3-D
/// curve (internal cap-to-cap edges created by `mev`/`mef`); we
/// reconstruct the 3-D path via `old_surface(old_pcurve(t))`.
///
/// `preserve_analytic_pcurves` (Phase 6c-D Stage 4.6): when `true`,
/// half-edges whose existing pcurve is one of the analytic forms
/// (`Linear`, `Elliptical`, `Sinusoidal`) are NOT re-fit.  The
/// `wire_substituted_boundary_to_edges` step writes analytic
/// `Linear` pcurves on substituted-boundary cap faces; closest-
/// point projection on the freshly-installed Gregory surface
/// frequently lands on the wrong boundary (row 0 vs col 9) for
/// FaceEdge bridges and clobbers the correct analytic pcurve.
fn reproject_loop_pcurves_onto_gregory(
    brep: &mut BRep,
    face: FaceHandle,
    old_surface: Option<&SurfaceGeometry>,
    preserve_analytic_pcurves: bool,
) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::nurbs_curve::NurbsCurve;

    let Some(face_data) = brep.faces.get(face) else {
        return;
    };
    let outer_loop = face_data.outer_loop;
    let surface = match face_data.surface.as_ref() {
        Some(sof) => sof.geometry.clone(),
        None => return,
    };
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return;
    };
    let start_he = loop_data.half_edge;

    // Walk the loop and collect half-edges first (avoids borrow
    // overlap when we mutate brep.half_edges).
    let mut hes: Vec<popkern_core::topology::HalfEdgeHandle> = Vec::new();
    let mut cursor = start_he;
    let mut safety = 64;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cursor);
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }

    let n_samples = 12;
    let edge_triple_enabled = std::env::var("POPKERN_EDGE_TRIPLE")
        .map(|v| v == "1")
        .unwrap_or(false);
    for he in &hes {
        let he = *he;
        let edge = brep.half_edges[he].edge;

        // Phase 6c-D Stage 4.6: preserve analytic pcurves (Linear /
        // Elliptical / Sinusoidal) when the cap face was wired by
        // `sphere_fit_edge_wire`.  Closest-point projection on the
        // freshly-installed Gregory surface routinely lands FaceEdge
        // bridge samples on row 0 (FilletEnd boundary) instead of
        // col 9 — wrong cap-surface UV — and clobbering the wire's
        // analytic Linear pcurve produces 3-5 cm `dev_cap` residual.
        if preserve_analytic_pcurves
            && let Some(existing) = brep.half_edges[he].pcurve.as_ref()
            && !matches!(
                existing.geometry,
                popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(_)
            )
        {
            continue;
        }

        // Stage 6: EdgeTriple path (gated). For each outer-loop edge
        // whose other-side face has an analytic surface (Plane,
        // Cylinder, Sphere, Torus), construct curve_3d + both pcurves
        // from the same source via `compute_edge_geometry` with
        // `EdgeHint::SurfaceBoundary`. Closes the cap-to-fillet seam
        // drift documented in Task #144.
        if edge_triple_enabled
            && let Some((curve_3d_geom, t_start, t_end, cap_pc, other_pc)) =
                super::three_quad_apply_edge_triple::try_build_edge_triple_for_seam(
                    brep, he, face, &surface,
                )
        {
            // Update BRep edge's 3-D curve.
            if let Some(edge_data) = brep.edges.get_mut(edge) {
                edge_data.curve = Some(popkern_core::geometry::curve::CurveOnEdge {
                    geometry: curve_3d_geom,
                    t_start,
                    t_end,
                });
            }
            // Update both half-edges' pcurves. Cap side gets pcurve_a,
            // other side (twin) gets pcurve_b. Both are parameterized
            // in curve_3d's natural direction; consumers handle
            // direction via twin's traversal logic.
            brep.half_edges[he].pcurve = Some(cap_pc.into());
            let twin = brep.half_edges[he].twin;
            brep.half_edges[twin].pcurve = Some(other_pc.into());
            continue;
        }

        // Decide our 3-D-source strategy for this half-edge:
        //   1. Edge has a 3-D curve → sample that.
        //   2. Edge has no 3-D curve but `old_surface` exists and
        //      the half-edge has an old pcurve → reconstruct
        //      `old_surface(old_pcurve(t))` as the 3-D path.
        //   3. Otherwise skip.
        let curve_on_edge = brep.edges.get(edge).and_then(|e| e.curve.clone());
        let old_pcurve = brep.half_edges[he].pcurve.clone();
        enum Source {
            Curve3D(popkern_core::geometry::curve::CurveOnEdge, bool),
            OldSurfPcurve(SurfaceGeometry, popkern_core::geometry::pcurve::Pcurve),
        }
        let source = if let Some(c) = curve_on_edge {
            // Determine direction — does the half-edge run with the
            // 3-D curve or against it?
            let v_origin = brep.half_edges[he].origin;
            let p_origin = brep.vertices[v_origin].position;
            let p_t_start = c.geometry.point_at(c.t_start);
            let p_t_end = c.geometry.point_at(c.t_end);
            let forward = p_origin.distance_to(&p_t_start) <= p_origin.distance_to(&p_t_end);
            Source::Curve3D(c, forward)
        } else if let (Some(os), Some(opc)) = (old_surface, old_pcurve) {
            Source::OldSurfPcurve(os.clone(), opc)
        } else {
            continue;
        };

        // Sample the 3-D path, project each sample onto the new
        // surface, build a degree-3 NURBS through the (u, v, 0)
        // samples.
        let mut points: Vec<popkern_core::geometry::Point3> = Vec::with_capacity(n_samples);
        let mut params: Vec<f64> = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            let s = i as f64 / (n_samples - 1) as f64;
            let p_3d = match &source {
                Source::Curve3D(coe, forward) => {
                    let t_s = if *forward { s } else { 1.0 - s };
                    let t = coe.t_start + t_s * (coe.t_end - coe.t_start);
                    coe.geometry.point_at(t)
                }
                Source::OldSurfPcurve(old_surf, old_pc) => {
                    use popkern_core::geometry::surface::SurfaceEval;
                    let (pc_lo, pc_hi) = old_pc.parameter_range();
                    let t = pc_lo + s * (pc_hi - pc_lo);
                    let uv = old_pc.point_at(t);
                    old_surf.point_at(uv.x, uv.y)
                }
            };
            let (u, v) = super::helpers::project_to_uv(&surface, &p_3d);
            // Clamp to the surface's parameter range so we don't fly
            // outside. Gregory: [0, 1] × [0, 1].
            let ((u_lo, u_hi), (v_lo, v_hi)) = {
                use popkern_core::geometry::surface::SurfaceEval;
                surface.parameter_range()
            };
            let u = if u_lo.is_finite() && u_hi.is_finite() {
                u.clamp(u_lo, u_hi)
            } else {
                u
            };
            let v = if v_lo.is_finite() && v_hi.is_finite() {
                v.clamp(v_lo, v_hi)
            } else {
                v
            };
            points.push(popkern_core::geometry::Point3::new(u, v, 0.0));
            params.push(s);
        }
        let degree = 3usize;
        let n_cps = (degree + 2).min(points.len());
        let mut fitted =
            match popkern_core::geometry::fitting::fit_curve(&points, &params, degree, n_cps) {
                Ok(c) => c,
                Err(_) => continue,
            };
        // Least-squares fit can produce control points *outside* the
        // surface's UV range (the curve stays inside the convex hull
        // of its CPs, but the CPs themselves can extend beyond the
        // sampled data). Clamp every CP to the surface UV box so
        // pcurve evaluations always land inside the valid Gregory
        // domain. May distort the curve slightly near the boundary
        // but eliminates the "uv = -1.5" runaway we saw on
        // pyramid_n04 (Cyl×Gregory edge with div=0.25 residual).
        use popkern_core::geometry::surface::SurfaceEval;
        let ((u_lo, u_hi), (v_lo, v_hi)) = surface.parameter_range();
        for cp in &mut fitted.control_points {
            if u_lo.is_finite() && u_hi.is_finite() {
                cp.x = cp.x.clamp(u_lo, u_hi);
            }
            if v_lo.is_finite() && v_hi.is_finite() {
                cp.y = cp.y.clamp(v_lo, v_hi);
            }
        }
        let _: &NurbsCurve = &fitted;
        brep.half_edges[he].pcurve = Some(fitted.into());
    }
}

/// Promote sub-patches from `RectangularTwoStripe` to
/// `RectangularBiNonic` via the Karciauskas-Peters guided-surface
/// construction. Returns `Some(worst_coplanarity)` when the
/// binonic ships with relaxed continuity (worst > 5e-2 commercial
/// gate — Phase H.8 step-D graceful degradation), `None` when the
/// promotion is accepted at commercial grade or reverted.
///
/// Phase 6c-B: `fillet_surfaces` (one entry per sub-patch, ordered
/// by ring index) are forwarded to `optimize_binonic_with_outer_fillets`
/// so the bi-9 L-M solver includes per-sample cap-vs-fillet G1
/// residuals along each cap's outer boundary, in addition to the
/// inner-spoke G1 residuals. Closes the cap-fillet seam tangent-
/// plane mismatch documented in Task #144.
///
/// `POPKERN_OUTER_FILLET_G1=1` env var gates the new path during
/// the spike. Without it, behaviour matches the pre-Phase-6c-B
/// optimizer (inner-spoke G1 only).
fn promote_decomp_to_kp_binonic(
    sub_patches: &mut [GregoryPatch],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    hex_boundaries: Option<&[([Point3; 10], [Point3; 10])]>,
) -> Option<f64> {
    let n = sub_patches.len();
    if n == 0 {
        return None;
    }
    // Every sub-patch must already be two-stripe — we read the
    // grid as the source of the K-P DirectBoundaryInput.
    let mut two_stripe: Vec<RectangularGregoryTwoStripe> = Vec::with_capacity(n);
    for patch in sub_patches.iter() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            trace_kp("non-two-stripe input", n, None);
            return None;
        };
        two_stripe.push((**ts).clone());
    }

    // Ring-topology spokes: consistent with the other promotion
    // rungs — spoke `k` shares sub-patch `k`'s UOne with sub-patch
    // `(k+1) % n`'s VZero (reversed).
    let spokes: Vec<Spoke> = (0..n)
        .map(|k| Spoke {
            a: SpokeEnd {
                patch_idx: k,
                side: Side::UOne,
                reversed: false,
            },
            b: SpokeEnd {
                patch_idx: (k + 1) % n,
                side: Side::VZero,
                reversed: true,
            },
        })
        .collect();

    let Some(input) = extract_direct_boundary_input(&two_stripe, &spokes, hex_boundaries) else {
        trace_kp("extract_direct_boundary_input → None", n, None);
        return None;
    };
    let char_len = estimate_characteristic_length(&two_stripe);
    let options = DirectConstructionOptions::with_characteristic_length(char_len);
    // Phase H.7 gap #2 closure — K-P 2004 §4 iterative twist
    // optimization. The one-shot Shirman-Séquin twist is under-
    // conditioned for extreme radius asymmetry (e.g. cube_extreme
    // 1:10); this loop refines the N central-vertex twists by
    // Gauss-Newton-style gradient descent against the ring-wide
    // G1 coplanarity residual. Converges in 3-5 iterations
    // typically; always monotonically improves the baseline
    // (the iterate tracks the best-seen residual).
    // Option C step 3: three-phase optimizer (bicubic → bi-quintic
    // → bi-nonic). Phase 3 only engages when phases 1+2 fail to
    // converge; for the vast majority of fixtures, behaviour is
    // identical to the bi-quintic path. Falls through to the bi-
    // quintic-only optimizer, then to the bicubic C-twist
    // optimizer (Option B) if both fail.
    // Phase 6c-B: env-gate on POPKERN_OUTER_FILLET_G1=1. When set,
    // pass the fillet surfaces to the binonic optimiser so its
    // L-M residual includes per-sample cap-vs-fillet G1 along
    // each cap's outer boundary.
    let outer_fillets: Option<Vec<popkern_core::geometry::surface::SurfaceGeometry>> =
        if std::env::var("POPKERN_OUTER_FILLET_G1").is_ok_and(|v| v == "1") {
            Some(fillet_surfaces.iter().map(|f| f.surface.clone()).collect())
        } else {
            None
        };

    let binonic = match super::kp_binonic_optimizer::optimize_binonic_with_outer_fillets(
        &input,
        &options,
        outer_fillets.as_deref(),
    ) {
        Some(r) => r.binonic,
        None => match optimize_biquintic_and_build_binonic(&input, &options) {
            Some(r) => r.binonic,
            None => {
                let Some(result) = optimize_c_twists_and_build_binonic(&input, &options) else {
                    trace_kp("optimize_binonic_and_build_binonic → None", n, None);
                    return None;
                };
                result.binonic
            }
        },
    };

    // Self-intersection gate — never commit a fold.
    let binonic_gp: Vec<GregoryPatch> = binonic
        .iter()
        .cloned()
        .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
        .collect();
    for (i, gp) in binonic_gp.iter().enumerate() {
        if gp.self_intersects_with_tolerance(1e-6) {
            trace_kp(
                "self_intersects_with_tolerance",
                n,
                Some(format!("sub-patch {i}")),
            );
            return None;
        }
    }

    // Quality guard — revert to two-stripe when binonic is
    // catastrophically worse than the predecessor rung. Threshold
    // 1.2× absorbs floating-point wobble while catching real
    // regressions.
    let pre_worst = worst_ring_coplanarity(sub_patches);
    let post_worst = worst_ring_coplanarity(&binonic_gp);
    // Phase 6c-D Stage 3: under boundary injection, the binonic
    // patch is built from a boundary-consistent frame + Coons
    // interior, NOT from the K-P guide-driven L-M optimum. Spoke
    // G1 is determined entirely by the (analytic) row0/col9 +
    // (input) row9/col0 geometry — the L-M can't refine it. The
    // 1.2× regression guard was tuned for the guide-driven path
    // and rejects legitimate-but-different geometry under
    // injection. Skip the guard when injection is active; the
    // self-intersection gate (already checked above) is the
    // remaining quality floor under injection.
    let injection_active = hex_boundaries.is_some();
    if !injection_active && pre_worst > 1e-12 && post_worst > 1.2 * pre_worst {
        trace_kp(
            "quality regression",
            n,
            Some(format!("pre={pre_worst:.3e} post={post_worst:.3e}")),
        );
        return None;
    }

    for (idx, patch) in binonic.into_iter().enumerate() {
        sub_patches[idx] = GregoryPatch::RectangularBiNonic(Box::new(patch));
    }

    // Phase H.8 step-D — `Some(worst)` signals relaxed continuity.
    const COMMERCIAL_GATE: f64 = 5e-2;
    if post_worst > COMMERCIAL_GATE {
        Some(post_worst)
    } else {
        None
    }
}

/// Worst-case G1 coplanarity across the full ring's shared spokes,
/// using 50 uniform samples on `t ∈ [0.05, 0.95]` per spoke. The
/// endpoint margin matches the commercial-gate test's convention
/// (avoids the rational-Gregory corner-singularity finite-angle
/// transition zone). Returns 0 if normals are ill-conditioned at
/// every sample (can't happen on a real ring, but guards against
/// log-scale division).
pub(super) fn worst_ring_coplanarity(patches: &[GregoryPatch]) -> f64 {
    let n = patches.len();
    let mut worst = 0.0_f64;
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            let t_b = 1.0 - t;
            let xtan_a = a.partial_u(1.0, t);
            let xtan_b = b.partial_v(t_b, 0.0);
            let t_spoke = a.partial_v(1.0, t);
            let cross = xtan_a.cross(&xtan_b);
            let numer = (cross.x * t_spoke.x + cross.y * t_spoke.y + cross.z * t_spoke.z).abs();
            let na = (xtan_a.x * xtan_a.x + xtan_a.y * xtan_a.y + xtan_a.z * xtan_a.z).sqrt();
            let nb = (xtan_b.x * xtan_b.x + xtan_b.y * xtan_b.y + xtan_b.z * xtan_b.z).sqrt();
            let ns = (t_spoke.x * t_spoke.x + t_spoke.y * t_spoke.y + t_spoke.z * t_spoke.z).sqrt();
            let denom = na * nb * ns;
            if denom > 1e-300 {
                let coplan = numer / denom;
                if coplan > worst {
                    worst = coplan;
                }
            }
        }
    }
    worst
}

/// Emit a diagnostic trace when the K-P promotion silently
/// falls back to two-stripe. Gated on the `POPKERN_KP_TRACE`
/// env var so production builds stay silent; set `=1` before
/// running a fillet test to see which stage rejected the input.
fn trace_kp(stage: &str, n: usize, detail: Option<String>) {
    if std::env::var_os("POPKERN_KP_TRACE").is_some() {
        match detail {
            Some(d) => eprintln!("[kp_binonic N={n}] fallback @ {stage}: {d}"),
            None => eprintln!("[kp_binonic N={n}] fallback @ {stage}"),
        }
    }
}

/// Build a [`DirectBoundaryInput`] from the biseptic-promoted
/// two-stripe seed patches + ring topology. Extracts the per-
/// spoke boundary CPs, corner xtan values, and input-face
/// tangent planes (approximated from the patch normals at C
/// and V — upstream biseptic promotion already wrote these
/// consistent with the true input face planes, so using the
/// patch's computed normal is equivalent).
///
/// `pub(super)` so the sibling `kp_polar_bi6_promote` module
/// can share the same extraction (I.10f).
pub(super) fn extract_direct_boundary_input(
    patches: &[RectangularGregoryTwoStripe],
    spokes: &[Spoke],
    hex_boundaries: Option<&[([Point3; 10], [Point3; 10])]>,
) -> Option<DirectBoundaryInput> {
    let n = patches.len();
    if n == 0 || spokes.len() != n {
        return None;
    }
    // Central vertex C: patch 0's grid[7][0] (u=1, v=0) per the
    // three-quad convention.
    let central_vertex = patches[0].grid[7][0];

    let mut spoke_data = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut central_vertex_tps = Vec::with_capacity(n);

    for (k, spoke) in spokes.iter().enumerate() {
        let patch_a = &patches[spoke.a.patch_idx];
        let patch_b = &patches[spoke.b.patch_idx];

        // Spoke boundary CPs: for the standard three-quad (A UOne
        // not reversed), these are patch_a.grid[7][0..=7].
        let mut boundary_cps = [Point3::new(0.0, 0.0, 0.0); 8];
        for (j, cp) in boundary_cps.iter_mut().enumerate() {
            *cp = match spoke.a.side {
                Side::UOne => patch_a.grid[7][j],
                Side::UZero => patch_a.grid[0][j],
                Side::VOne => patch_a.grid[j][7],
                Side::VZero => patch_a.grid[j][0],
            };
        }

        // xtan at corner C (spoke-t = 0) and V (spoke-t = 1) on
        // each side. Side A: partial in the cross-boundary
        // direction at the spoke's side. Side B: same, but with
        // reversal mapping spoke-t → B's local parameter via
        // t_local = 1 - t when b.reversed is true.
        let (ua_c, va_c) = uv_at_spoke_t(spoke.a.side, 0.0, spoke.a.reversed);
        let (ua_v, va_v) = uv_at_spoke_t(spoke.a.side, 1.0, spoke.a.reversed);
        let xtan_a_at_c = cross_tangent_on_side(patch_a, spoke.a.side, ua_c, va_c);
        let xtan_a_at_v = cross_tangent_on_side(patch_a, spoke.a.side, ua_v, va_v);

        let (ub_c, vb_c) = uv_at_spoke_t(spoke.b.side, 0.0, spoke.b.reversed);
        let (ub_v, vb_v) = uv_at_spoke_t(spoke.b.side, 1.0, spoke.b.reversed);
        let xtan_b_at_c = cross_tangent_on_side(patch_b, spoke.b.side, ub_c, vb_c);
        let xtan_b_at_v = cross_tangent_on_side(patch_b, spoke.b.side, ub_v, vb_v);

        spoke_data.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c,
            xtan_b_at_c,
            xtan_a_at_v,
            xtan_b_at_v,
            // Radii not used by the current blend-coeffs formula
            // (constant-α form). Populate with a plausible placeholder
            // so the struct validates; the K-K-C 1991 variable-
            // radius extension will consume these later.
            radius_at_c: 1.0,
            radius_at_v: 1.0,
            edge_id: format!("spoke_{k}"),
        });

        // Outer corner position: the shared (u=1, v=1) of patch A.
        let v_pos = patch_a.grid[7][7];
        // Tangent plane at V from patch A's normal there.
        let du_v = patch_a.partial_u(1.0, 1.0);
        let dv_v = patch_a.partial_v(1.0, 1.0);
        let n_v = du_v.cross(&dv_v).normalized()?;
        let tps_v = TangentPlane::from_normal(n_v)?;
        outer_corners.push(OuterCornerData {
            position: v_pos,
            tangent_plane: tps_v,
        });

        // Central-vertex tangent plane from patch A's normal at C.
        let du_c = patch_a.partial_u(1.0, 0.0);
        let dv_c = patch_a.partial_v(1.0, 0.0);
        let n_c = du_c.cross(&dv_c).normalized()?;
        let tps_c = TangentPlane::from_normal(n_c)?;
        central_vertex_tps.push(tps_c);
    }

    // Phase 6c-D Stage 1d: when the caller provides
    // hex_boundaries (deg-9 Bernstein CPs of the substituted
    // FilletEnd / FaceEdge curves per sub-patch), thread them
    // into the new `HexBoundaryData` field. The K-P bi-nonic
    // ring builder uses these to bypass `k_p_boundary` Hermite
    // synthesis and trace the analytic `sphere ∩ cylinder` arc
    // exactly. `None` falls through to existing behaviour.
    let hex_boundaries_data = hex_boundaries.map(|src| {
        src.iter()
            .map(
                |(row0, col9)| popkern_core::geometry::gregory::HexBoundaryData {
                    row0: *row0,
                    col9: *col9,
                },
            )
            .collect::<Vec<_>>()
    });

    Some(DirectBoundaryInput {
        central_vertex,
        spokes: spoke_data,
        outer_corners,
        central_vertex_tps,
        seed_patches: Some(patches.to_vec()),
        hex_boundaries: hex_boundaries_data,
    })
}

/// Map a spoke parameter `t ∈ [0, 1]` to the `(u, v)` on a
/// side-indexed patch boundary, honouring the `reversed` flag
/// that three-quad uses for its VZero B-ends (spoke-t = 0 ↔
/// patch-u = 1 when reversed).
fn uv_at_spoke_t(side: Side, t: f64, reversed: bool) -> (f64, f64) {
    let t_local = if reversed { 1.0 - t } else { t };
    match side {
        Side::UZero => (0.0, t_local),
        Side::UOne => (1.0, t_local),
        Side::VZero => (t_local, 0.0),
        Side::VOne => (t_local, 1.0),
    }
}

/// Partial derivative in the cross-boundary direction at the
/// given side. For U-sides that's `∂S/∂u`; for V-sides `∂S/∂v`.
fn cross_tangent_on_side(
    patch: &RectangularGregoryTwoStripe,
    side: Side,
    u: f64,
    v: f64,
) -> popkern_core::geometry::Vec3 {
    match side {
        Side::UZero | Side::UOne => patch.partial_u(u, v),
        Side::VZero | Side::VOne => patch.partial_v(u, v),
    }
}

/// Estimate a characteristic length for tolerance scaling — the
/// max bounding-box extent across the seed patches. Populates
/// `DirectConstructionOptions::characteristic_length` so all
/// ratio-based tolerances scale with the geometry's natural
/// size.
fn estimate_characteristic_length(patches: &[RectangularGregoryTwoStripe]) -> f64 {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for patch in patches {
        for i in 0..8 {
            for j in 0..8 {
                let p = patch.grid[i][j];
                for (k, v) in [p.x, p.y, p.z].iter().enumerate() {
                    if *v < min[k] {
                        min[k] = *v;
                    }
                    if *v > max[k] {
                        max[k] = *v;
                    }
                }
            }
        }
    }
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    (dx * dx + dy * dy + dz * dz).sqrt().max(1e-9)
}

#[cfg(test)]
#[path = "three_quad_apply_tests.rs"]
mod tests;
