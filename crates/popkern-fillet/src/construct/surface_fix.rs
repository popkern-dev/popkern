//! Surface correction utilities for multi-edge fillets.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

/// Correct face surface assignments after topology modifications.
///
/// Checks every face. If a face has no surface, or if its surface does
/// not match its boundary vertices (max vertex-to-surface distance
/// exceeds threshold), replaces it with the best matching surface from
/// the pool of original planar surfaces and fillet surfaces.
pub(crate) fn correct_face_surfaces(
    brep: &mut BRep,
    original: &BRep,
    fillet_surfaces: &[SurfaceGeometry],
) {
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceOnFace};
    use popkern_core::geometry::uv_projection::project_point_to_uv;

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

        // Check if current surface is geometrically correct.
        let current_max_dist = brep.faces[fh]
            .surface
            .as_ref()
            .map(|sof| {
                positions
                    .iter()
                    .filter_map(|p| {
                        project_point_to_uv(&sof.geometry, p)
                            .map(|(u, v)| p.distance_to(&sof.geometry.point_at(u, v)))
                    })
                    .fold(0.0_f64, f64::max)
            })
            .unwrap_or(f64::INFINITY);

        if current_max_dist < 0.01 {
            continue;
        }

        // Find the best matching surface that improves on current.
        // Never replace a non-planar surface with a plane (fillet
        // surfaces are intentionally assigned even if vertices aren't
        // exactly on them).
        let current_is_nonplanar = brep.faces[fh]
            .surface
            .as_ref()
            .is_some_and(|s| !matches!(s.geometry, SurfaceGeometry::Plane { .. }));

        let mut best: Option<(SurfaceOnFace, f64)> = None;
        for sof in &candidates {
            // Don't downgrade non-planar to planar.
            if current_is_nonplanar && matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
                continue;
            }
            let max_dist: f64 = positions
                .iter()
                .map(|p| {
                    project_point_to_uv(&sof.geometry, p)
                        .map(|(u, v)| p.distance_to(&sof.geometry.point_at(u, v)))
                        .unwrap_or(f64::INFINITY)
                })
                .fold(0.0_f64, f64::max);

            if max_dist < current_max_dist {
                match &best {
                    Some((_, bd)) if max_dist >= *bd => {}
                    _ => best = Some((sof.clone(), max_dist)),
                }
            }
        }

        // Try a vertex blend sphere: if the face is small (≤6 verts)
        // and its centroid is near a corner of the original solid,
        // a sphere at (corner ± R) with radius R might fit.
        if best.as_ref().is_none_or(|(_, d)| *d > 0.08)
            && positions.len() <= 6
            && let Some(sphere_sof) = try_sphere_for_corner(&positions, original, fillet_surfaces)
        {
            let sd: f64 = positions
                .iter()
                .map(|p| match &sphere_sof.geometry {
                    SurfaceGeometry::Sphere { center, radius } => {
                        let dist = ((p.x - center.x).powi(2)
                            + (p.y - center.y).powi(2)
                            + (p.z - center.z).powi(2))
                        .sqrt();
                        (dist - radius).abs()
                    }
                    _ => f64::INFINITY,
                })
                .fold(0.0_f64, f64::max);
            if sd < current_max_dist {
                match &best {
                    Some((_, bd)) if sd >= *bd => {}
                    _ => best = Some((sphere_sof, sd)),
                }
            }
        }

        if let Some((sof, dist)) = best {
            // Only assign if the candidate actually fits reasonably well.
            // Without this threshold, a wrong surface that happens to be
            // the "least bad" in the pool would still be assigned.
            if dist < 0.5 {
                brep.faces[fh].surface = Some(sof);
                continue;
            }
        }

        // Last resort: fit a plane through face vertices. Only for faces
        // that don't already have a non-planar surface (fillet surfaces
        // are intentionally assigned even if imprecise).
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

    // Average fillet radius from torus minor radii.
    let avg_radius = fillet_surfaces
        .iter()
        .filter_map(|s| match s {
            SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
            _ => None,
        })
        .next()
        .unwrap_or(0.1);

    // Face centroid.
    let n = positions.len() as f64;
    let cx = positions.iter().map(|p| p.x).sum::<f64>() / n;
    let cy = positions.iter().map(|p| p.y).sum::<f64>() / n;
    let cz = positions.iter().map(|p| p.z).sum::<f64>() / n;

    // Find nearest original vertex.
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
