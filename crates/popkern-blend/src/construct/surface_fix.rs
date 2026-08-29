//! Surface correction utilities for multi-edge fillets.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, TopologyRole};

use super::strip_orientation::reconcile_if_reauthored_blend_surface;

/// Compute the *max* distance from a list of boundary positions to
/// a candidate surface. Used as the "is this surface acceptable at
/// all" gate, before falling through to candidate replacement.
fn max_fit_score(geom: &SurfaceGeometry, positions: &[popkern_core::geometry::Point3]) -> f64 {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    positions
        .iter()
        .filter_map(|p| {
            project_point_to_uv(geom, p).map(|(u, v)| p.distance_to(&geom.point_at(u, v)))
        })
        .fold(0.0_f64, f64::max)
}

/// Correct face surface assignments after topology modifications.
///
/// Checks every face. If a face has no surface, or if its surface does
/// not match its boundary vertices (median vertex-to-surface distance
/// exceeds threshold), replaces it with the best matching surface from
/// the pool of original planar surfaces and fillet surfaces.
pub(crate) fn correct_face_surfaces(
    brep: &mut BRep,
    original: &BRep,
    fillet_surfaces: &[SurfaceGeometry],
) {
    use popkern_core::geometry::surface::SurfaceOnFace;

    let mut candidates: Vec<SurfaceOnFace> = original
        .faces
        .iter()
        .filter_map(|(_, f)| f.surface.clone())
        .collect();
    for surf in fillet_surfaces {
        candidates.push(SurfaceOnFace {
            geometry: surf.clone(),
            orientation: true,
        });
    }

    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();

    for fh in face_handles {
        if brep.faces[fh].role == TopologyRole::ConstructionResidual {
            continue;
        }
        let positions = brep.face_outer_loop_positions(fh);
        if positions.is_empty() {
            continue;
        }

        // Check if the current surface is acceptable as-is. Use
        // **max** distance for the skip gate (so a face whose
        // surface fits *all* its vertices well is preserved
        // unchanged) AND for candidate comparison (so we don't
        // accidentally pick a surface that perfectly fits the
        // majority but has huge errors at outlier vertices). The
        // unconditional vertex-blend sphere fallback below is the
        // escape hatch for cases where the cube-plane "wins" max
        // fit but the semantically correct surface is the sphere.
        let current_max = brep.faces[fh]
            .surface
            .as_ref()
            .map(|sof| max_fit_score(&sof.geometry, &positions))
            .unwrap_or(f64::INFINITY);

        if current_max < 0.01 {
            continue;
        }

        // the collapsed-projection fallback work / the sphere-fit edge-wiring work: trust upstream Sphere and Gregory
        // assignments. Per-fillet Sphere sub-caps (cube_mixed/extreme,
        // pyramid_uniform with rung ON) and Gregory sub-caps
        // (rung OFF, the sphere-fit edge-wiring work path) are assigned with intent in
        // `three_quad_apply`. Their boundaries can include inner-spoke
        // vertices that lie on *other* sub-face surfaces, so
        // `current_max` against this sub-face's surface can be ~r —
        // larger than the 1cm gate above. The cube-corner-inscribed-
        // sphere fallback in `try_sphere_for_corner` would clobber
        // both with the wrong cube-style centre.
        // NURBS caps are deliberately NOT trusted here: a cube-style
        // corner's NURBS blend is legitimately replaced by the exact
        // rolling-ball sphere; try_sphere_for_corner's orthogonality
        // guard (PK-221) protects non-orthogonal polar caps instead.
        let current_is_cap_surface = brep.faces[fh].surface.as_ref().is_some_and(|s| {
            matches!(
                s.geometry,
                SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Gregory(_)
            )
        });
        if current_is_cap_surface {
            continue;
        }

        // Trust CONSTRUCTED FILLET STRIP surfaces the same way:
        // assigned with intent; corner ANCHOR vertices ~r·(√2−1) off
        // the strip (pre-route-B corner topology) fail the 1cm gate,
        // and the refit swapped real vertical strips for wall planes
        // — the corner passes own those vertices, not this surface.
        let current_is_fillet_surface = brep.faces[fh].surface.as_ref().is_some_and(|s| {
            fillet_surfaces.iter().any(|fs| {
                super::sliver_absorb_classifiers::surfaces_equal(
                    &s.geometry,
                    fs,
                    &popkern_core::tolerance::Tolerance::default(),
                )
            })
        });
        if current_is_fillet_surface {
            continue;
        }

        let current_score = current_max;

        // Find the best matching surface that improves on current.
        // Never replace a non-planar surface with a plane (fillet
        // surfaces are intentionally assigned even if vertices aren't
        // exactly on them).
        let current_is_nonplanar = brep.faces[fh]
            .surface
            .as_ref()
            .is_some_and(|s| !matches!(s.geometry, SurfaceGeometry::Plane { .. }));

        let is_sphere_eligible = positions.len() <= 6;
        let mut best: Option<(SurfaceOnFace, f64)> = None;
        for sof in &candidates {
            // Don't downgrade non-planar to planar.
            if current_is_nonplanar && matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
                continue;
            }
            // For sphere-eligible faces, don't let fillet surfaces
            // (Cylinder/Torus) compete — defer to the sphere check.
            if is_sphere_eligible
                && matches!(
                    sof.geometry,
                    SurfaceGeometry::Cylinder { .. } | SurfaceGeometry::Torus { .. }
                )
            {
                continue;
            }
            let score = max_fit_score(&sof.geometry, &positions);
            if score < current_score {
                match &best {
                    Some((_, bd)) if score >= *bd => {}
                    _ => best = Some((sof.clone(), score)),
                }
            }
        }

        // Try a vertex blend sphere unconditionally for small
        // faces. A vertex-blend face's boundary often has the
        // 3-4 sphere contact corners + 1-3 cube-corner outliers,
        // and the cube corners *also* lie on cube planes — so a
        // cube plane wins the median fit even though the
        // semantically correct surface is the blend sphere. The
        // threshold-gated fallback above missed these because the
        // plane's median was low. We now try the sphere for every
        // small face (≤6 verts) and let `try_sphere_for_corner`'s
        // own "centroid near a cube corner" filter prune real
        // plane faces.
        //
        // When the small Plane face is adjacent
        // to a Sphere or Gregory cap, skip the sphere fallback.
        // Class D sliver_absorb is the right path — it merges the
        // Plane fragment INTO the cap. Without this guard,
        // pyramid_uniform's apex Plane corner fragments get
        // prematurely converted to cube-style inscribed Spheres
        // (geometrically wrong for pyramid) before sliver_absorb
        // can run.
        let face_adjacent_to_cap_surface = is_adjacent_to_sphere_or_gregory(brep, fh);
        if positions.len() <= 6
            && !face_adjacent_to_cap_surface
            && let Some(sphere_sof) = try_sphere_for_corner(&positions, original, fillet_surfaces)
        {
            let sd = max_fit_score(&sphere_sof.geometry, &positions);
            // SCALE-RELATIVE acceptance: the pre-repair cap loop
            // contains the stale cube-corner outlier at distance
            // exactly r(√2−1) from the true sphere (the apex vertex
            // relocation repair runs AFTER surface assignment), so
            // the fit score grows linearly with the fillet radius
            // and a 0.5-absolute gate silently stops promoting the
            // exact sphere above r ≈ 1.2 (measured: uniform r=3
            // cube corner — sd = 1.2426 = r(√2−1), the Gregory
            // NURBS shipped instead, and its boundary lied 0.0712·r
            // off the strips; check-9 catch). The earlier 0.5 floor
            // is kept so sub-unit radii keep their historical
            // margin.
            let radius_scale = fillet_surfaces
                .iter()
                .filter_map(crate::blend::sphere_fit_centre::fillet_nominal_radius)
                .fold(1.0_f64, f64::max);
            let sphere_gate = 0.5 * radius_scale;
            if popkern_core::trace::trace_enabled("sphere-fit") {
                eprintln!(
                    "[surface-fix] face {fh:?} sphere candidate sd={sd:.4e} gate={sphere_gate:.4e}"
                );
            }
            // Prefer a sphere for small corner faces whenever one
            // fits cleanly. The exact-analytic sphere beats any
            // NURBS approximation (including the §5.3 Gregory
            // polar patch) for representation clarity, even when
            // the NURBS passes through every boundary sample by
            // construction. The `try_sphere_for_corner` helper's
            // own centroid-near-cube-corner filter prunes faces
            // where a sphere is semantically wrong.
            if sd < sphere_gate {
                let prefer_sphere = match &best {
                    Some((existing, bd)) => {
                        let existing_is_non_sphere =
                            !matches!(existing.geometry, SurfaceGeometry::Sphere { .. });
                        sd < *bd || (existing_is_non_sphere && sd <= *bd * 1.5)
                    }
                    None => true,
                };
                if prefer_sphere {
                    best = Some((sphere_sof, sd));
                }
            }
        }

        if let Some((sof, dist)) = best {
            // Only assign if the candidate actually fits reasonably
            // well (same radius-scaled window as the sphere gate —
            // the stale-outlier distance grows with r for every
            // candidate kind).
            let assign_gate = 0.5
                * fillet_surfaces
                    .iter()
                    .filter_map(crate::blend::sphere_fit_centre::fillet_nominal_radius)
                    .fold(1.0_f64, f64::max);
            if dist < assign_gate {
                brep.faces[fh].surface = Some(sof);
                reconcile_if_reauthored_blend_surface(brep, fh, fillet_surfaces);
                continue;
            }
        }

        // Last resort: fit a plane through face vertices. Only for
        // faces that don't already have a non-planar surface
        // (fillet surfaces are intentionally assigned even if
        // imprecise).
        if current_is_nonplanar {
            continue;
        }
        if positions.len() >= 3 {
            let p0 = positions[0];
            let p1 = positions[1];
            let p2 = positions[2];
            let u = popkern_core::geometry::Vec3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
            let v = popkern_core::geometry::Vec3::new(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
            let n = u.cross(&v);
            let len = n.length();
            if len > 1e-12 {
                let normal = popkern_core::geometry::Vec3::new(n.x / len, n.y / len, n.z / len);
                let ul = u.length();
                let u_axis = popkern_core::geometry::Vec3::new(u.x / ul, u.y / ul, u.z / ul);
                let v_axis = normal.cross(&u_axis);
                brep.faces[fh].surface = Some(SurfaceOnFace {
                    geometry: SurfaceGeometry::Plane {
                        origin: p0,
                        normal,
                        u_axis,
                        v_axis,
                    },
                    orientation: true,
                });
            }
        }
    }
}

/// Returns `true` if any half-edge in `face`'s outer
/// loop has its twin face's surface set to Sphere or Gregory.
/// Used as a guard against the small-Plane → Sphere fallback in
/// `correct_face_surfaces`: when a Plane corner fragment is
/// adjacent to an upstream cap surface, Class D sliver_absorb
/// merges it into the cap, which is the canonical path. The
/// `try_sphere_for_corner` fallback would prematurely convert it
/// to a (geometrically wrong) Sphere on pyramid_uniform.
fn is_adjacent_to_sphere_or_gregory(brep: &BRep, face: FaceHandle) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return false;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut safety = 16usize;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        if let Some(twin_loop) = brep.half_edges.get(he.twin).map(|t| t.loop_)
            && let Some(twin_face) = brep.loops.get(twin_loop).map(|l| l.face)
            && twin_face != face
            && let Some(twin_sof) = brep.faces.get(twin_face).and_then(|f| f.surface.as_ref())
            && matches!(
                twin_sof.geometry,
                SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Gregory(_)
            )
        {
            return true;
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    false
}

/// Try to construct a sphere surface for a vertex blend corner.
///
/// Finds the nearest vertex of the original solid. The sphere center
/// is at (vertex ± R) in the direction of each adjacent face's inward
/// normal, with radius R = average fillet radius.
fn try_sphere_for_corner(
    positions: &[popkern_core::geometry::Point3],
    original: &BRep,
    fillet_surfaces: &[SurfaceGeometry],
) -> Option<popkern_core::geometry::surface::SurfaceOnFace> {
    use popkern_core::geometry::surface::SurfaceOnFace;

    // Fillet radii from torus minor / cylinder radii. For the
    // sphere to be semantically correct as a corner patch, all
    // incident fillet strips must share one radius — otherwise
    // the rolling ball isn't the same at each contact and a
    // Gregory patch (or analogous N-sided construction) is the
    // right answer. Bail out on mixed radii so `build_smooth_blend`
    // keeps whatever N-sided NURBS it produced.
    // A swept NURBS strip is not a constant circular section — it is
    // a variable-radius or conic-profile blend (PK-189), and the
    // rolling ball that would justify a corner sphere does not exist
    // for it. Promoting anyway shipped a Sphere cap whose boundary
    // arcs lay the conic shoulder gap off the strips (measured
    // 2.7e-3 CurveSurfaceDeviation on the rho=0.35 cube corner).
    if fillet_surfaces
        .iter()
        .any(|s| matches!(s, SurfaceGeometry::Nurbs(_)))
    {
        return None;
    }
    let radii: Vec<f64> = fillet_surfaces
        .iter()
        .filter_map(|s| match s {
            SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
            SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    if radii.len() >= 2 {
        let (min_r, max_r) = radii.iter().fold((f64::INFINITY, 0.0_f64), |(lo, hi), &r| {
            (lo.min(r), hi.max(r))
        });
        if max_r - min_r > 1e-6 {
            return None;
        }
    }
    let avg_radius = radii.first().copied().unwrap_or(0.1);

    // Orthogonality guard (PK-221): the centre construction below
    // offsets PER GLOBAL AXIS from the corner, which is only
    // meaningful for a cube-style corner — mutually orthogonal,
    // axis-aligned fillet axes. A pyramid apex (non-orthogonal
    // supports) satisfied every other gate and shipped a sphere with
    // g1 ~1.0 at its rim. Decline so the constructed N-sided cap
    // stands.
    let axes: Vec<popkern_core::geometry::Vec3> = fillet_surfaces
        .iter()
        .filter_map(|s| match s {
            SurfaceGeometry::Torus { axis, .. } | SurfaceGeometry::Cylinder { axis, .. } => {
                axis.normalized()
            }
            _ => None,
        })
        .collect();
    const AXIS_ALIGNMENT_TOLERANCE: f64 = 1.0e-6;
    for axis in &axes {
        let max_component = axis.x.abs().max(axis.y.abs()).max(axis.z.abs());
        if max_component < 1.0 - AXIS_ALIGNMENT_TOLERANCE {
            return None;
        }
    }
    for (i, a) in axes.iter().enumerate() {
        for b in axes.iter().skip(i + 1) {
            if a.dot(b).abs() > AXIS_ALIGNMENT_TOLERANCE
                && a.dot(b).abs() < 1.0 - AXIS_ALIGNMENT_TOLERANCE
            {
                return None;
            }
        }
    }

    // Use the per-coordinate **median** position rather than the
    // mean centroid. The boundary of a multi-fillet vertex blend
    // can include outlier vertices that bleed in from adjacent
    // fillets' regions; the mean centroid is biased by them, but
    // the median sits inside the dominant cluster and is robust.
    let median_axis = |xs: &mut Vec<f64>| -> f64 {
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        xs[xs.len() / 2]
    };
    let mut xs: Vec<f64> = positions.iter().map(|p| p.x).collect();
    let mut ys: Vec<f64> = positions.iter().map(|p| p.y).collect();
    let mut zs: Vec<f64> = positions.iter().map(|p| p.z).collect();
    let cx = median_axis(&mut xs);
    let cy = median_axis(&mut ys);
    let cz = median_axis(&mut zs);

    // Find nearest original vertex to the median position.
    let mut best_dist = f64::INFINITY;
    let mut best_pos = None;
    for (_, v) in &original.vertices {
        let d = ((v.position.x - cx).powi(2)
            + (v.position.y - cy).powi(2)
            + (v.position.z - cz).powi(2))
        .sqrt();
        if d < best_dist {
            best_dist = d;
            best_pos = Some(v.position);
        }
    }

    let corner = best_pos?;
    if best_dist > avg_radius * 3.0 {
        return None; // Too far from any corner.
    }

    // Sphere center: offset from corner toward the solid interior.
    // For each axis, move by R in the direction away from the nearest
    // face (toward the interior).
    //
    // **Direction-check via centroid, not median.** The median sits
    // inside the dominant cluster (robust to outliers) but on a
    // symmetric N-vertex cap (e.g. cube_uniform's 3 rolling-ball
    // footprints at distance r along each cube axis), the median
    // equals the corner position exactly — defeating the
    // `>`-comparison and forcing the wrong (`else`) branch, which
    // produced sphere centres OUTSIDE the solid (cube_uniform
    // shipped with centre = (-r, -r, -r) instead of (+r, +r, +r),
    // explaining the long-standing 5-15 cm cap-vs-cube drift).
    //
    // The centroid (mean) breaks the symmetric tie consistently:
    // for the same 3 footprints, centroid = (r/3, r/3, r/3),
    // strictly greater than the corner on every axis. We use the
    // median for `cx/cy/cz` (it's the right choice for *position*),
    // but the centroid for the *direction* of the offset.
    let mean_x = positions.iter().map(|p| p.x).sum::<f64>() / positions.len() as f64;
    let mean_y = positions.iter().map(|p| p.y).sum::<f64>() / positions.len() as f64;
    let mean_z = positions.iter().map(|p| p.z).sum::<f64>() / positions.len() as f64;
    let sx = if mean_x > corner.x {
        corner.x + avg_radius
    } else {
        corner.x - avg_radius
    };
    let sy = if mean_y > corner.y {
        corner.y + avg_radius
    } else {
        corner.y - avg_radius
    };
    let sz = if mean_z > corner.z {
        corner.z + avg_radius
    } else {
        corner.z - avg_radius
    };

    Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Sphere {
            center: popkern_core::geometry::Point3::new(sx, sy, sz),
            radius: avg_radius,
        },
        orientation: true,
    })
}
