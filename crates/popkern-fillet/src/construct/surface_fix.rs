//! Surface correction utilities for multi-edge fillets.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

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
        if positions.len() <= 6
            && let Some(sphere_sof) = try_sphere_for_corner(&positions, original, fillet_surfaces)
        {
            let sd = max_fit_score(&sphere_sof.geometry, &positions);
            if sd < current_score {
                let prefer_sphere = match &best {
                    Some((existing, bd)) => {
                        let existing_is_non_sphere =
                            !matches!(existing.geometry, SurfaceGeometry::Sphere { .. });
                        // Prefer the sphere over planes and fillet
                        // surfaces (Cylinder/Torus) when the sphere
                        // score is within 50% of the candidate score.
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
            // Only assign if the candidate actually fits reasonably well.
            if dist < 0.5 {
                brep.faces[fh].surface = Some(sof);
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

    // Average fillet radius from torus minor radii or cylinder radii.
    // Straight-edge plane-plane fillets produce Cylinder surfaces;
    // curved-edge fillets produce Torus surfaces.
    let avg_radius = fillet_surfaces
        .iter()
        .filter_map(|s| match s {
            SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
            SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
            _ => None,
        })
        .next()
        .unwrap_or(0.1);

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
    let sx = if cx > corner.x {
        corner.x + avg_radius
    } else {
        corner.x - avg_radius
    };
    let sy = if cy > corner.y {
        corner.y + avg_radius
    } else {
        corner.y - avg_radius
    };
    let sz = if cz > corner.z {
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
