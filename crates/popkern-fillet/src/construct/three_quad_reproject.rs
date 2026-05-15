//! Gregory pcurve reprojection for newly assigned cap sub-faces.

use super::*;
use crate::construct::{helpers, three_quad_apply_edge_triple};
use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::topology::{SurfaceBoundaryHintDecline, SurfaceBoundaryHintDeclineReason};

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
pub(super) fn reproject_loop_pcurves_onto_gregory(
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

        // Stage 6: EdgeTriple path (gated). For each outer-loop edge
        // whose other-side face has an analytic surface (Plane,
        // Cylinder, Sphere, Torus), construct curve_3d + both pcurves
        // from the same source via `compute_edge_geometry` with
        // `EdgeHint::SurfaceBoundary`. Closes the cap-to-fillet seam
        // drift documented in Task #144.
        if edge_triple_enabled
            && let Some((curve_3d_geom, t_start, t_end, cap_pc, other_pc)) =
                three_quad_apply_edge_triple::try_build_edge_triple_for_seam(
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
            let cap_pcurve = cap_pc.into();
            record_boundary_evidence_from_pcurve(brep, he, &cap_pcurve);
            brep.half_edges[he].pcurve = Some(cap_pcurve);
            let twin = brep.half_edges[he].twin;
            brep.half_edges[twin].pcurve = Some(other_pc.into());
            continue;
        }

        // Phase 6c-D Stage 4.6: preserve analytic pcurves (Linear /
        // Elliptical / Sinusoidal) when the cap face was wired by
        // `sphere_fit_edge_wire`.  Closest-point projection on the
        // freshly-installed Gregory surface routinely lands FaceEdge
        // bridge samples on row 0 (FilletEnd boundary) instead of
        // col 9 — wrong cap-surface UV — and clobbering the wire's
        // analytic Linear pcurve produces 3-5 cm `dev_cap` residual.
        if preserve_analytic_pcurves
            && let Some(existing) = brep.half_edges[he].pcurve.clone()
            && !matches!(
                existing.geometry,
                popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(_)
            )
        {
            if matches!(surface, SurfaceGeometry::Gregory(_)) {
                record_boundary_evidence_from_pcurve(brep, he, &existing);
            }
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
            let (u, v) = helpers::project_to_uv(&surface, &p_3d);
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
        let fitted_pcurve = popkern_core::geometry::pcurve::Pcurve::from(fitted);
        if matches!(surface, SurfaceGeometry::Gregory(_)) {
            record_boundary_evidence_from_pcurve(brep, he, &fitted_pcurve);
        }
        brep.half_edges[he].pcurve = Some(fitted_pcurve);
    }
    if matches!(surface, SurfaceGeometry::Gregory(_)) {
        record_unclassified_gregory_outer_loop_boundaries(brep, &hes);
    }
}

fn record_boundary_evidence_from_pcurve(
    brep: &mut BRep,
    half_edge: popkern_core::topology::HalfEdgeHandle,
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
) {
    let classification = classify_rectangular_boundary_pcurve(pcurve);
    if let Some(boundary) = classification.boundary {
        brep.set_surface_boundary_hint(half_edge, boundary);
    } else {
        brep.set_surface_boundary_decline(
            half_edge,
            SurfaceBoundaryHintDecline {
                reason: SurfaceBoundaryHintDeclineReason::NotOnCanonicalBoundary,
                best_boundary_deviation: classification.best_deviation,
                best_boundary: classification.boundary,
                endpoint_gap: 0.0,
            },
        );
    }
}

struct BoundaryPcurveClassification {
    boundary: Option<BoundarySide>,
    best_deviation: f64,
}

fn classify_rectangular_boundary_pcurve(
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
) -> BoundaryPcurveClassification {
    let (t0, t1) = pcurve.parameter_range();
    let mut samples = Vec::with_capacity(9);
    for i in 0..9 {
        let s = i as f64 / 8.0;
        samples.push(pcurve.point_at(t0 + s * (t1 - t0)));
    }
    classify_rectangular_boundary_samples(&samples)
}

fn classify_rectangular_boundary_samples(
    samples: &[popkern_core::geometry::Point3],
) -> BoundaryPcurveClassification {
    let threshold = 1.0e-6;
    let candidates = [
        (
            BoundarySide::UMin,
            samples
                .iter()
                .map(|point| point.x.abs())
                .fold(0.0, f64::max),
        ),
        (
            BoundarySide::UMax,
            samples
                .iter()
                .map(|point| (1.0 - point.x).abs())
                .fold(0.0, f64::max),
        ),
        (
            BoundarySide::VMin,
            samples
                .iter()
                .map(|point| point.y.abs())
                .fold(0.0, f64::max),
        ),
        (
            BoundarySide::VMax,
            samples
                .iter()
                .map(|point| (1.0 - point.y).abs())
                .fold(0.0, f64::max),
        ),
    ];
    let Some(best) = candidates.into_iter().min_by(|a, b| a.1.total_cmp(&b.1)) else {
        return BoundaryPcurveClassification {
            boundary: None,
            best_deviation: f64::INFINITY,
        };
    };
    BoundaryPcurveClassification {
        boundary: (best.1 <= threshold).then_some(best.0),
        best_deviation: best.1,
    }
}

fn record_unclassified_gregory_outer_loop_boundaries(
    brep: &mut BRep,
    half_edges: &[popkern_core::topology::HalfEdgeHandle],
) {
    for &half_edge in half_edges {
        if brep.surface_boundary_hint(half_edge).is_some()
            || brep.surface_boundary_decline(half_edge).is_some()
        {
            continue;
        }
        let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(half_edge_data.twin) else {
            continue;
        };
        let Some(twin_face) = brep.loops.get(twin.loop_).map(|loop_data| loop_data.face) else {
            continue;
        };
        let analytic_neighbor = brep
            .faces
            .get(twin_face)
            .and_then(|face| face.surface.as_ref())
            .is_some_and(|surface| {
                matches!(
                    surface.geometry,
                    SurfaceGeometry::Plane { .. }
                        | SurfaceGeometry::Cylinder { .. }
                        | SurfaceGeometry::Sphere { .. }
                        | SurfaceGeometry::Torus { .. }
                )
            });
        if !analytic_neighbor {
            continue;
        }
        let Some(pcurve) = brep.half_edges[half_edge].pcurve.clone() else {
            brep.set_surface_boundary_decline(
                half_edge,
                SurfaceBoundaryHintDecline {
                    reason: SurfaceBoundaryHintDeclineReason::MissingPcurve,
                    best_boundary_deviation: f64::INFINITY,
                    best_boundary: None,
                    endpoint_gap: 0.0,
                },
            );
            continue;
        };
        record_boundary_evidence_from_pcurve(brep, half_edge, &pcurve);
    }
}

/// Verify that every half-edge on a cap face's outer trim loop owns a
/// pcurve after surface assignment and reprojection.
pub(super) fn verify_outer_loop_pcurves(brep: &BRep, face: FaceHandle) -> Result<(), String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| "face handle is invalid".to_string())?;
    let loop_data = brep
        .loops
        .get(face_data.outer_loop)
        .ok_or_else(|| "outer loop handle is invalid".to_string())?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..brep.half_edges.len().max(1) {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| "outer loop half-edge handle is invalid".to_string())?;
        if he.pcurve.is_none() {
            return Err(format!("half-edge {cursor:?} has no pcurve"));
        }
        cursor = he.next;
        if cursor == start {
            return Ok(());
        }
    }
    Err("outer loop did not close while checking pcurves".to_string())
}

/// Fill missing outer-loop pcurves from the live half-edge endpoint
/// vertices after a cap face has received its final surface.
///
/// This is intentionally narrow: it only authors pcurves that are
/// absent, uses the face's current surface, and leaves existing
/// analytic or fitted pcurves untouched. The caller still runs
/// `verify_outer_loop_pcurves` and the surrounding transaction rolls
/// back if the face remains invalid.
pub(super) fn author_missing_outer_loop_pcurves_from_vertices(brep: &mut BRep, face: FaceHandle) {
    let Some(face_data) = brep.faces.get(face) else {
        return;
    };
    let outer_loop = face_data.outer_loop;
    let Some(surface) = face_data.surface.as_ref().map(|sof| sof.geometry.clone()) else {
        return;
    };
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..brep.half_edges.len().max(1) {
        let Some(he) = brep.half_edges.get(cursor).cloned() else {
            return;
        };
        if he.pcurve.is_none()
            && let Some(next) = brep.half_edges.get(he.next)
            && let (Some(origin), Some(end)) =
                (brep.vertices.get(he.origin), brep.vertices.get(next.origin))
        {
            let start_uv = helpers::project_to_uv(&surface, &origin.position);
            let end_uv = helpers::project_to_uv(&surface, &end.position);
            brep.half_edges[cursor].pcurve = Some(popkern_core::geometry::pcurve::Pcurve::linear(
                start_uv, end_uv, 0.0, 1.0,
            ));
        }
        cursor = he.next;
        if cursor == start {
            return;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn verify_outer_loop_pcurves_accepts_complete_face() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        assert!(verify_outer_loop_pcurves(&brep, face).is_ok());
    }

    #[test]
    fn verify_outer_loop_pcurves_rejects_missing_trim_pcurve() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let he = brep.loops[brep.faces[face].outer_loop].half_edge;
        brep.half_edges[he].pcurve = None;

        let err = verify_outer_loop_pcurves(&brep, face).unwrap_err();

        assert!(err.contains("has no pcurve"), "{err}");
    }

    #[test]
    fn author_missing_outer_loop_pcurves_fills_only_missing_trims() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let he = brep.loops[brep.faces[face].outer_loop].half_edge;
        let next = brep.half_edges[he].next;
        let existing_next_was_some = brep.half_edges[next].pcurve.is_some();
        brep.half_edges[he].pcurve = None;

        author_missing_outer_loop_pcurves_from_vertices(&mut brep, face);

        assert!(brep.half_edges[he].pcurve.is_some());
        assert_eq!(
            brep.half_edges[next].pcurve.is_some(),
            existing_next_was_some
        );
        assert!(verify_outer_loop_pcurves(&brep, face).is_ok());
    }
}
