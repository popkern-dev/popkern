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

/// The "already fits" window as a fraction of the smallest fillet
/// radius (see `correct_face_surfaces`).
const FIT_WINDOW_PER_RADIUS: f64 = 0.05;
/// The historical absolute window, kept as the cap so unit-scale
/// fixtures (r >= 0.2) see no change.
const FIT_WINDOW_MAX: f64 = 0.01;

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

    // "Already fits" window for a face's current surface. A corner's
    // pre-repair cap loop carries the stale cube-corner outlier r(√2−1)
    // = 0.414 r off the true surface, so the window must sit well under
    // that at every scale. It was an absolute 0.01: below r ≈ 0.024 the
    // outlier hid inside it, the corner cap kept its wrong surface and
    // the fillet declined with every corner edge 0.414 r off its faces
    // (a 0.7 mm round on a 16 mm cube, i.e. any mm-scale part; cubes
    // from 0.008 to 0.25 failed, 0.3 and up served). 0.05 r keeps the
    // historical 0.01 wherever r >= 0.2.
    let smallest_radius = fillet_surfaces
        .iter()
        .filter_map(crate::blend::sphere_fit_centre::fillet_nominal_radius)
        .fold(f64::INFINITY, f64::min);
    let fit_window = (smallest_radius * FIT_WINDOW_PER_RADIUS).min(FIT_WINDOW_MAX);

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

        if current_max < fit_window {
            continue;
        }

        // the collapsed-projection fallback work / the sphere-fit edge-wiring work: trust upstream Sphere and Gregory
        // assignments. Per-fillet Sphere sub-caps (cube_mixed/extreme,
        // pyramid_uniform with rung ON) and Gregory sub-caps
        // (rung OFF, the sphere-fit edge-wiring work path) are assigned with intent in
        // `three_quad_apply`. Their boundaries can include inner-spoke
        // vertices that lie on *other* sub-face surfaces, so
        // `current_max` against this sub-face's surface can be ~r —
        // larger than the fit window above. The cube-corner-inscribed-
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
        // the strip (pre-route-B corner topology) fail the fit window,
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

#[path = "surface_fix_sphere.rs"]
mod sphere;
use sphere::try_sphere_for_corner;
