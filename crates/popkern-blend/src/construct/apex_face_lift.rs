//! Lift apex-region pyramid face surfaces to
//! contain the substituted FaceEdge bridge.
//!
//! After sphere-fit boundary substitution, the FaceEdge bridge
//! between adjacent FilletEnd terminations sits 2-4 cm off the
//! original pyramid face plane (sphere∩cyl arc terminations are
//! at z ≈ 0.815 while the rolling-ball footprint sat at z ≈
//! 0.755).  The cap surface contains the bridge exactly (Stage
//! 4.6); the original `Plane` face cannot — every plane UV maps
//! to a point on the plane, never to the off-plane bridge.
//!
//! This module replaces the affected pyramid face's `Plane`
//! surface with a degree-(1, 3) Bezier patch:
//!   * `u = 0` row: 4 control points along the substituted bridge
//!   * `u = 1` row: 4 control points on the original plane, at the
//!     orthogonal projection of the bridge CPs
//!
//! The surface bilinearly interpolates from off-plane (u=0,
//! containing the bridge) to on-plane (u=1).  Other edges of the
//! pyramid face (base, spines) are away from the apex and project
//! into u ≈ 1 region of the new surface, near the original plane.
//!
//! Always on. Closes `cap_outer_g0` on `pyramid4_mixed` from ~3 cm
//! to near zero on the bridge edges.

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Replace `face`'s `Plane` surface with a degree-(1, 3) Bezier
/// patch whose `u=0` boundary traces the substituted bridge curve
/// on edge `bridge_edge`.  Returns `Ok(())` on success; `Err(...)`
/// when the face has no plane surface, the bridge edge has no
/// curve, or some other precondition fails.
///
/// After the surface swap, callers must refit the face's outer-
/// loop pcurves onto the new NURBS surface (the existing pcurves
/// are in the old plane's UV space and become invalid).
pub(super) fn lift_apex_face_with_bridge(
    brep: &mut BRep,
    face: FaceHandle,
    bridge_edge: EdgeHandle,
) -> Result<(), String> {
    use popkern_core::geometry::curve::CurveEval;

    // Pull the original support surface. Plane supports use an exact
    // orthogonal projection; already-lifted or repaired NURBS supports use
    // closest-point projection into the previous surface.
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| "face not found".to_string())?;
    let sof = face_data
        .surface
        .as_ref()
        .ok_or_else(|| "face has no surface".to_string())?;
    let source_surface = sof.geometry.clone();
    let plane_projection = if let SurfaceGeometry::Plane { origin, normal, .. } = &source_surface {
        Some((
            *origin,
            normal
                .normalized()
                .ok_or_else(|| "plane normal is zero".to_string())?,
        ))
    } else {
        None
    };

    // Pull the bridge's deg-9 Bernstein samples (10 points at
    // uniform parameters) and downsample to 4 CPs by picking
    // t = 0, 1/3, 2/3, 1.  The bridge is a straight line so a
    // degree-3 Bezier through these is exact.
    let edge_data = brep
        .edges
        .get(bridge_edge)
        .ok_or_else(|| "bridge edge not found".to_string())?;
    let coe = edge_data
        .curve
        .as_ref()
        .ok_or_else(|| "bridge edge has no curve3d".to_string())?;
    let (t0, t1) = (coe.t_start, coe.t_end);
    let span = t1 - t0;
    if span.abs() < 1e-12 {
        return Err("bridge edge has zero parameter span".into());
    }
    let bridge_cps: [Point3; 4] = [
        coe.geometry.point_at(t0),
        coe.geometry.point_at(t0 + span / 3.0),
        coe.geometry.point_at(t0 + 2.0 * span / 3.0),
        coe.geometry.point_at(t1),
    ];

    // Project each bridge CP onto the previous support surface; these become
    // the `u=1` row of the new patch.
    let project_to_plane = |p: Point3| -> Point3 {
        if let Some((plane_origin, plane_normal)) = plane_projection {
            let d = p - plane_origin;
            let h = d.dot(&plane_normal);
            return Point3::new(
                p.x - plane_normal.x * h,
                p.y - plane_normal.y * h,
                p.z - plane_normal.z * h,
            );
        }
        source_surface
            .closest_point(&p)
            .map(|(u, v)| source_surface.point_at(u, v))
            .unwrap_or_else(|| source_surface.point_at(0.0, 0.0))
    };
    let plane_cps: [Point3; 4] = [
        project_to_plane(bridge_cps[0]),
        project_to_plane(bridge_cps[1]),
        project_to_plane(bridge_cps[2]),
        project_to_plane(bridge_cps[3]),
    ];

    // Build a degree-(1, 3) Bezier surface.  Row-major layout:
    // `cps[i_u * num_v + i_v]` with num_u = 2, num_v = 4.
    //   row u=0: bridge_cps        (off-plane, i_u=0)
    //   row u=1: plane_cps         (on-plane,  i_u=1)
    let mut control_points = Vec::with_capacity(8);
    control_points.extend_from_slice(&bridge_cps);
    control_points.extend_from_slice(&plane_cps);
    let weights = vec![1.0; 8];
    let knots_u = vec![0.0, 0.0, 1.0, 1.0];
    let knots_v = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];

    let surface = NurbsSurface {
        degree_u: 1,
        degree_v: 3,
        knots_u,
        knots_v,
        control_points,
        num_u: 2,
        num_v: 4,
        weights,
    };

    // Swap the face's surface.
    let mut new_sof = sof.clone();
    new_sof.geometry = SurfaceGeometry::Nurbs(surface);
    brep.faces[face].surface = Some(new_sof);

    // Project each bridge CP back gives a point on the plane;
    // double-check the construction is non-degenerate.
    let max_lift: f64 = bridge_cps
        .iter()
        .zip(plane_cps.iter())
        .map(|(b, p)| {
            let d = *b - *p;
            (d.x * d.x + d.y * d.y + d.z * d.z).sqrt()
        })
        .fold(0.0_f64, f64::max);
    if max_lift < 1e-9 {
        // Bridge already on plane — no real change, but the swap
        // is still applied (NURBS encoding of the plane).  Caller
        // can choose to skip the lift entirely.
    }

    Ok(())
}

/// Walk every cap face's outer loop and lift each adjacent Plane
/// face that's not already lifted.  Refits the lifted face's
/// outer-loop pcurves onto the new NURBS surface so subsequent
/// pipeline stages see consistent UV mappings.  Returns the
/// number of plane faces lifted.
pub(crate) fn lift_all_apex_plane_faces(brep: &mut BRep, cap_faces: &[FaceHandle]) -> usize {
    let trace = popkern_core::trace::trace_enabled("apex-face-lift");
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let mut already_lifted: std::collections::BTreeSet<FaceHandle> =
        std::collections::BTreeSet::new();
    // Discover (plane_face, bridge_edge) pairs by walking each cap face's outer loop.
    let mut jobs: Vec<(FaceHandle, EdgeHandle)> = Vec::new();
    for &cap_face in cap_faces {
        let Some(face_data) = brep.faces.get(cap_face) else {
            continue;
        };
        let outer_loop = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start_he = loop_data.half_edge;
        let mut cur = start_he;
        let mut step = 0usize;
        #[allow(clippy::while_let_loop)]
        loop {
            let Some(he) = brep.half_edges.get(cur) else {
                break;
            };
            let twin_face = brep
                .loops
                .get(brep.half_edges[he.twin].loop_)
                .map(|l| l.face);
            if let Some(tf) = twin_face
                && !cap_set.contains(&tf)
                && let Some(twin_sof) = brep.faces.get(tf).and_then(|f| f.surface.as_ref())
                && matches!(
                    twin_sof.geometry,
                    SurfaceGeometry::Plane { .. }
                        | SurfaceGeometry::Cylinder { .. }
                        | SurfaceGeometry::Nurbs(_)
                )
                && !already_lifted.contains(&tf)
            {
                if trace {
                    eprintln!(
                        "[apex-face-lift] cap={cap_face:?} step={step} plane={tf:?} bridge_edge={:?}",
                        he.edge
                    );
                }
                jobs.push((tf, he.edge));
                already_lifted.insert(tf);
            }
            cur = he.next;
            step += 1;
            if cur == start_he || step > 32 {
                break;
            }
        }
    }
    if trace {
        eprintln!("[apex-face-lift] {} jobs queued", jobs.len());
    }
    // Apply each lift, then refit the lifted face's outer-loop pcurves.
    let mut count = 0;
    for (plane_face, bridge_edge) in jobs {
        match lift_apex_face_with_bridge(brep, plane_face, bridge_edge) {
            Ok(()) => {
                if trace {
                    eprintln!("[apex-face-lift]   OK plane={plane_face:?} bridge={bridge_edge:?}");
                }
                refit_face_outer_loop_pcurves(brep, plane_face);
                // After surface swap, the bridge edge's plane-side
                // pcurve must be analytic Linear (0, 0) → (0, 1) —
                // the bridge IS the v-axis boundary of the lifted
                // patch by construction.  `refit_face_outer_loop_pcurves`
                // uses Gauss-Newton, which on this degenerate-Jacobian
                // patch (`partial_u` is constant = plane-minus-bridge,
                // `partial_v` varies) converges to corner (1, 1) or
                // (1, 0) instead of (0, s).  Direct analytic pcurve
                // bypasses the issue.
                write_analytic_bridge_pcurve_on_plane_side(brep, plane_face, bridge_edge);
                count += 1;
            }
            Err(msg) => {
                if trace {
                    eprintln!(
                        "[apex-face-lift]   ERR plane={plane_face:?} bridge={bridge_edge:?} : {msg}"
                    );
                }
            }
        }
    }
    count
}

pub(super) fn write_analytic_boundary_pcurve_on_support_side(
    brep: &mut BRep,
    support_face: FaceHandle,
    bridge_edge: EdgeHandle,
    u_value: f64,
) {
    use popkern_core::geometry::pcurve::Pcurve;
    let Some(edge_data) = brep.edges.get(bridge_edge) else {
        return;
    };
    let Some(coe) = edge_data.curve.as_ref() else {
        return;
    };
    let (t_start, t_end) = (coe.t_start, coe.t_end);
    let primary = edge_data.half_edge;
    let twin = brep.half_edges[primary].twin;
    let primary_face = brep.loops[brep.half_edges[primary].loop_].face;
    let support_he = if primary_face == support_face {
        primary
    } else {
        twin
    };
    let support_face_check = brep.loops[brep.half_edges[support_he].loop_].face;
    if support_face_check != support_face {
        return;
    }

    let v_origin = brep.half_edges[support_he].origin;
    let p_origin = brep.vertices[v_origin].position;
    let p_t_start = {
        use popkern_core::geometry::curve::CurveEval;
        coe.geometry.point_at(t_start)
    };
    let p_t_end = {
        use popkern_core::geometry::curve::CurveEval;
        coe.geometry.point_at(t_end)
    };
    let forward = p_origin.distance_to(&p_t_start) <= p_origin.distance_to(&p_t_end);
    let (start_uv, end_uv) = if forward {
        ((u_value, 0.0), (u_value, 1.0))
    } else {
        ((u_value, 1.0), (u_value, 0.0))
    };
    let pc = Pcurve::linear(start_uv, end_uv, t_start, t_end);
    brep.half_edges[support_he].pcurve = Some(pc);
}

/// Write an analytic Linear pcurve `(0, 0) → (0, 1)` on the
/// plane-side half-edge of `bridge_edge`.  After the lift,
/// `surface(0, v) == bridge.point_at(t_start + v · span)` exactly,
/// so the analytic pcurve gives bit-exact correspondence between
/// the cap-side and plane-side surface evaluations.
fn write_analytic_bridge_pcurve_on_plane_side(
    brep: &mut BRep,
    plane_face: FaceHandle,
    bridge_edge: EdgeHandle,
) {
    write_analytic_boundary_pcurve_on_support_side(brep, plane_face, bridge_edge, 0.0);
}

/// Re-fit every outer-loop half-edge's pcurve on `face` against
/// the face's current surface.  Used after the surface swap to
/// bring stale (in-the-old-Plane-UV-space) pcurves into the new
/// NURBS UV space.
pub(super) fn refit_face_outer_loop_pcurves(brep: &mut BRep, face: FaceHandle) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::fitting::fit_curve;
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    let Some(face_data) = brep.faces.get(face) else {
        return;
    };
    let surface = match face_data.surface.as_ref() {
        Some(sof) => sof.geometry.clone(),
        None => return,
    };
    let outer_loop = face_data.outer_loop;
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return;
    };
    let start_he = loop_data.half_edge;
    let mut hes: Vec<popkern_core::topology::HalfEdgeHandle> = Vec::new();
    let mut cur = start_he;
    let mut safety = 64usize;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cur);
        let Some(he) = brep.half_edges.get(cur) else {
            break;
        };
        cur = he.next;
        if cur == start_he {
            break;
        }
    }
    let n_samples = 12;
    let ((u_lo, u_hi), (v_lo, v_hi)) = surface.parameter_range();
    for he in hes {
        let edge_h = brep.half_edges[he].edge;
        let Some(coe) = brep.edges.get(edge_h).and_then(|e| e.curve.clone()) else {
            // No 3-D curve → leave existing pcurve alone.
            continue;
        };
        // Direction sense: does this half-edge run with the curve?
        let v_origin = brep.half_edges[he].origin;
        let p_origin = brep.vertices[v_origin].position;
        let p_t_start = coe.geometry.point_at(coe.t_start);
        let p_t_end = coe.geometry.point_at(coe.t_end);
        let forward = p_origin.distance_to(&p_t_start) <= p_origin.distance_to(&p_t_end);
        let mut points = Vec::with_capacity(n_samples);
        let mut params = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            let s = i as f64 / (n_samples - 1) as f64;
            let t_s = if forward { s } else { 1.0 - s };
            let t = coe.t_start + t_s * (coe.t_end - coe.t_start);
            let p_3d = coe.geometry.point_at(t);
            let Some((u, v)) = project_point_to_uv(&surface, &p_3d) else {
                continue;
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
            points.push(Point3::new(u, v, 0.0));
            params.push(s);
        }
        if points.len() < 2 {
            continue;
        }
        let degree = 3usize.min(points.len() - 1);
        let n_cps = (degree + 2).min(points.len());
        let Ok(fitted) = fit_curve(&points, &params, degree, n_cps) else {
            continue;
        };
        brep.half_edges[he].pcurve = Some(fitted.into());
    }
}
