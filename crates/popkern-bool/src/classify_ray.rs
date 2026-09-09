//! Single-ray transversal crossing count for the classify stage's
//! point-in-solid voting (split from `classify.rs` under the
//! 500-line file budget). The edge-band poisoning, cross-face
//! t-dedup, and extension-band rules live here.

use super::*;

#[path = "classify_ray_gates.rs"]
mod gates;
use gates::HitGate;

// ---------------------------------------------------------------------------
// 3. Single-ray crossing count
// ---------------------------------------------------------------------------

/// Cast a single ray and count transversal crossings with the faces of a solid.
///
/// Returns `Ok(count)` with the number of transversal face crossings, or
/// `Err(())` if any intersection is ambiguous (near-tangent or near an edge).
pub(super) fn ray_crossing_count(
    origin: &Point3,
    direction: &Vec3,
    brep: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
    skip_hole_hits: bool,
) -> Result<usize, ()> {
    let ray_curve = CurveGeometry::Line {
        origin: *origin,
        direction: *direction,
    };
    let t_range = (0.0, 1e6);
    // Accepted crossing parameters, deduplicated ACROSS faces before
    // the parity count: a crossing through (or near) a shared edge is
    // one crossing of the closed body but is legitimately reported by
    // BOTH adjacent faces — summing per-face counts read it twice and
    // flipped the parity (measured on a thread coil: an inside-void
    // sample's clean +z ray counted its single flank crossing once
    // per adjacent wall, voted even, and the opening fragment was
    // kept — PK-460).
    let mut crossing_ts: Vec<f64> = Vec::new();

    let ray_trace = popkern_core::trace::trace_enabled("ray-classify");
    for (face_handle, aabb) in face_aabbs {
        if !ray_intersects_aabb(origin, direction, aabb) {
            if ray_trace {
                eprintln!("  [ray] {face_handle:?} aabb miss {aabb:?}");
            }
            continue;
        }

        let face = &brep.faces[*face_handle];
        let surface_on_face = match face.surface.as_ref() {
            Some(s) => s,
            None => continue,
        };
        let surface = &surface_on_face.geometry;

        let ((u_min, u_max), (v_min, v_max)) = face_parameter_range(brep, *face_handle);

        // NURBS faces query a PADDED window and treat band hits as an
        // inconclusive RAY: a ray grazing along a tool's corner edge
        // converges just outside both adjacent charts (measured on a
        // thread coil: the exit crossings landed at u = −0.014 and
        // u = 1.04 of the two walls, both range-rejected, the parity
        // read 0, and an inside-tool sample classified Outside —
        // PK-460). A crossing that near a chart edge cannot be
        // honestly assigned to either face; the vote loop's other
        // directions decide instead. Analytic faces keep their tuned
        // trim gates (below) untouched.
        let is_nurbs = matches!(
            surface,
            popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
        );
        let band = 0.05;
        let (qu, qv) = if is_nurbs {
            let (du, dv) = (u_max - u_min, v_max - v_min);
            (
                (u_min - band * du, u_max + band * du),
                (v_min - band * dv, v_max + band * dv),
            )
        } else {
            ((u_min, u_max), (v_min, v_max))
        };
        let result = intersect_surface_curve(surface, qu, qv, &ray_curve, t_range, tolerance)
            .map_err(|_| ())?;
        if ray_trace {
            eprintln!(
                "  [ray] {face_handle:?} u({u_min:.2},{u_max:.2}) v({v_min:.2},{v_max:.2}) hits {}",
                result.intersections.len()
            );
        }

        for isect in &result.intersections {
            match isect {
                SurfaceCurveIntersection::Point { t, u, v, .. } => {
                    // Only forward intersections (t > 0).
                    if *t < tolerance.linear() {
                        continue;
                    }

                    // Chart-edge band verdict for NURBS faces (see the
                    // padded query above): a hit within the edge band
                    // — EITHER side of the native window — poisons
                    // the ray. One-sided rejection was measured
                    // nondeterministic on a thread tool's end-cap
                    // corner: the same ray's Newton landed at
                    // v = 1.0012 in one call (rejected) and
                    // v = 1.0000005 in the next (counted), flipping
                    // the parity per call. Within the band the
                    // crossing cannot be honestly assigned; a
                    // different direction decides.
                    if is_nurbs {
                        let (du, dv) = (u_max - u_min, v_max - v_min);
                        let band_u = 0.005 * du;
                        let band_v = 0.005 * dv;
                        let near_edge = (*u - u_min).abs() < band_u
                            || (*u - u_max).abs() < band_u
                            || (*v - v_min).abs() < band_v
                            || (*v - v_max).abs() < band_v;
                        let in_native = *u >= u_min + band_u
                            && *u <= u_max - band_u
                            && *v >= v_min + band_v
                            && *v <= v_max - band_v;
                        if near_edge || !in_native {
                            if in_native {
                            } else if *u < u_min - band_u
                                || *u > u_max + band_u
                                || *v < v_min - band_v
                                || *v > v_max + band_v
                            {
                                // Deep in the pad, far outside the
                                // face: an extended-basis phantom,
                                // not an edge ambiguity — skip it
                                // without poisoning the ray.
                                continue;
                            }
                            if popkern_core::trace::trace_enabled("ray-classify") {
                                eprintln!(
                                    "  [ray-band] {face_handle:?} u={u:.4} v={v:.4} in edge band — ray inconclusive"
                                );
                            }
                            return Err(());
                        }
                    }

                    let normal = surface.normal_at(*u, *v);
                    let n_len = normal.length();
                    if n_len < tolerance.numerical_floor() {
                        // Degenerate normal — skip.
                        continue;
                    }
                    let n_unit = normal * (1.0 / n_len);
                    let d_len = direction.length();
                    if d_len < tolerance.numerical_floor() {
                        return Err(());
                    }
                    let d_unit = *direction * (1.0 / d_len);

                    let dot = n_unit.dot(&d_unit).abs();

                    // If the hit is nearly tangent, the ray is inconclusive.
                    if dot < tolerance.angular() {
                        return Err(());
                    }

                    // A hit inside one of the face's INNER LOOPS is
                    // a hole, not a crossing: the parameter box is the
                    // trim proxy and covers holes (measured: a point
                    // in a counterbore void classified Inside because
                    // the ray "hit" the holed top face through its
                    // own hole; the tool's proud band was kept and the
                    // result went non-manifold). Scoped to PLANE
                    // faces with inner loops — the general trim-aware
                    // gate was measured to break simple cube CSG, and
                    // holes on periodic charts have unresolved seam
                    // branches.
                    let is_plane = matches!(
                        surface,
                        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
                    );
                    // Boundary-derived geometric gates (classify_ray_gates.rs):
                    // cylinder angular span, sphere cap, cylinder/cone axial
                    // span (FINAL bodies only), and the plane in-plane box
                    // (both modes). Evaluated in the original order with the
                    // original short-circuits: a SkipHit stops the chain (the
                    // hit is discarded before later gates run), PoisonRay
                    // aborts the ray. Each gate returns Accept immediately
                    // for non-matching surface types.
                    let verdict = if skip_hole_hits {
                        let mut v = gates::cylinder_angular_gate(
                            brep,
                            *face_handle,
                            surface,
                            origin,
                            direction,
                            *t,
                        );
                        if matches!(v, HitGate::Accept) {
                            v = gates::sphere_cap_gate(
                                brep,
                                *face_handle,
                                surface,
                                origin,
                                direction,
                                *t,
                                tolerance,
                            );
                        }
                        if matches!(v, HitGate::Accept) {
                            v = gates::axial_span_gate(
                                brep,
                                *face_handle,
                                surface,
                                origin,
                                direction,
                                *t,
                                tolerance,
                            );
                        }
                        v
                    } else {
                        HitGate::Accept
                    };
                    let verdict = if matches!(verdict, HitGate::Accept) {
                        gates::plane_box_gate(
                            brep,
                            *face_handle,
                            surface,
                            origin,
                            direction,
                            *t,
                            tolerance,
                        )
                    } else {
                        verdict
                    };
                    match verdict {
                        HitGate::Accept => {}
                        HitGate::SkipHit => continue,
                        HitGate::PoisonRay => return Err(()),
                    }
                    // TORUS hits: accept only within the face's chart
                    // CELL (PK-304). The torus's finite native range
                    // made the window the WHOLE surface, so one
                    // geometric crossing counted once per patch face —
                    // parity always even, every fragment classified
                    // Outside, and the union retained both operands
                    // whole (measured: 22 of 22 fragments kept,
                    // validate_all flagged 7 retained interior
                    // interpenetrations). Fold-containment against the
                    // cell is the trim truth; a hit within a sliver of
                    // the cell border is ambiguous (it belongs to the
                    // neighbouring patch as much as this one) and
                    // poisons the ray like a near-tangent hit.
                    if matches!(
                        surface,
                        popkern_core::geometry::surface::SurfaceGeometry::Torus { .. }
                    ) {
                        match crate::clip::torus_hit_in_face_cell(
                            brep,
                            *face_handle,
                            *u,
                            *v,
                            tolerance,
                        ) {
                            Some(true) => {}
                            Some(false) => {
                                if ray_trace {
                                    eprintln!(
                                        "  [ray-rej-cell] {face_handle:?} uv=({u:.3},{v:.3})"
                                    );
                                }
                                continue;
                            }
                            None => return Err(()),
                        }
                    }
                    // A hit OUTSIDE the face's outer pcurve polygon is
                    // not on the face at all — the padded parameter box
                    // is a solver window, not a trim (measured: a ray
                    // passing 6 mm UNDER a pocket wall crossed its
                    // padded box and a buried tool fragment classified
                    // Outside, dropping two pocket walls). Scoped to
                    // PLANE charts (affine pcurves, no seam branches)
                    // and to FINAL bodies (`skip_hole_hits`): mid-split
                    // states carry STALE pcurves — recompute runs at
                    // the end of the operation — and a stale polygon
                    // rejects hits that ARE on the face (measured:
                    // cube ∩ cube misclassified during the stage
                    // classification).
                    if skip_hole_hits
                        && is_plane
                        && !hit_inside_outer_loop(brep, *face_handle, *u, *v)
                    {
                        if popkern_core::trace::trace_enabled("ray-classify") {
                            eprintln!("  [ray-rej] face {face_handle:?} uv=({u:.2},{v:.2})");
                        }
                        continue;
                    }
                    if skip_hole_hits
                        && is_plane
                        && hit_inside_inner_loop(brep, *face_handle, *u, *v)
                    {
                        continue;
                    }
                    if popkern_core::trace::trace_enabled("ray-classify") {
                        let h = *origin + *direction * *t;
                        eprintln!(
                            "  [ray-acc] {face_handle:?} t={t:.2} hit=({:.2},{:.2},{:.2})",
                            h.x, h.y, h.z
                        );
                    }
                    // Transversal crossing.
                    crossing_ts.push(*t);
                }
                SurfaceCurveIntersection::Overlap { .. } => {
                    // Ray lies on the surface — inconclusive.
                    return Err(());
                }
            }
        }
    }

    // Distinct-crossing count: hits within the shared-edge slop of an
    // already-counted parameter are the same crossing seen from the
    // other face. Direction magnitudes are ~1, so parameter distance
    // is world distance.
    crossing_ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let dedup = tolerance.linear() * 10.0;
    let mut count = 0usize;
    let mut last: Option<f64> = None;
    for t in crossing_ts {
        if last.is_none_or(|l| t - l > dedup) {
            count += 1;
            last = Some(t);
        }
    }

    Ok(count)
}
