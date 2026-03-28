//! Surface correction and trim clipping utilities for multi-edge fillets.

use std::collections::HashMap;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

pub(crate) fn resolve<H: Copy + Eq + std::hash::Hash>(remap: &HashMap<H, H>, mut h: H) -> H {
    for _ in 0..64 {
        match remap.get(&h) {
            Some(&new) => h = new,
            None => return h,
        }
    }
    h
}

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

        if let Some((sof, _dist)) = best {
            brep.faces[fh].surface = Some(sof);
            continue;
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

/// Clip trim anchors to the face boundary. If an anchor is NOT on the
/// face (consumed by a previous fillet), find where the trim line crosses
/// a boundary edge, split it, and return the intersection vertex.
pub(crate) fn clip_anchors_to_face(
    brep: &mut BRep,
    face: FaceHandle,
    anchors: (VertexHandle, VertexHandle),
    tolerance: &Tolerance,
) -> (VertexHandle, VertexHandle) {
    let a = clip_one_anchor(brep, face, anchors.0, anchors.1, tolerance);
    let b = clip_one_anchor(brep, face, anchors.1, anchors.0, tolerance);
    (a, b)
}

/// Clip one anchor to the face boundary if it's not on the loop.
fn clip_one_anchor(
    brep: &mut BRep,
    face: FaceHandle,
    anchor: VertexHandle,
    other: VertexHandle,
    _tolerance: &Tolerance,
) -> VertexHandle {
    if is_vertex_on_face(brep, face, anchor) {
        return anchor;
    }

    let anchor_pos = brep.vertices[anchor].position;
    let other_pos = brep.vertices[other].position;

    // Walk the face boundary. For each edge, compute the intersection
    // of the trim line (anchor→other) with the edge line segment.
    let face_data = match brep.faces.get(face) {
        Some(f) => f,
        None => return anchor,
    };
    let start_he = match brep.loops.get(face_data.outer_loop) {
        Some(l) => l.half_edge,
        None => return anchor,
    };

    // Strategy 1: Find where the trim line crosses a boundary edge.
    let mut best_crossing: Option<(VertexHandle, Option<(EdgeHandle, f64)>)> = None;
    let mut best_cross_dist = f64::INFINITY;
    let mut cursor = start_he;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = &brep.half_edges[cursor];
        let next_he = &brep.half_edges[he.next];
        let he_origin = he.origin;
        let he_dest = next_he.origin;
        let e0 = brep.vertices[he_origin].position;
        let e1 = brep.vertices[he_dest].position;

        if let Some((t_edge, pt)) = line_line_nearest(anchor_pos, other_pos, e0, e1) {
            let dist = pt.distance_to(&anchor_pos);
            if dist < best_cross_dist {
                if t_edge <= 0.001 {
                    best_crossing = Some((he_origin, None));
                } else if t_edge >= 0.999 {
                    best_crossing = Some((he_dest, None));
                } else {
                    best_crossing = Some((he_origin, Some((he.edge, t_edge))));
                }
                best_cross_dist = dist;
            }
        }

        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }

    if let Some((vertex, split_info)) = best_crossing {
        return match split_info {
            None => vertex,
            Some((edge, t)) => brep.split_edge(edge, t).map_or(anchor, |s| s.vertex),
        };
    }

    // No crossing found — keep the original anchor.
    // The face split will fail, which is handled as a per-edge failure.
    anchor
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

/// Check if a vertex is on a face's outer loop boundary.
fn is_vertex_on_face(brep: &BRep, face: FaceHandle, vertex: VertexHandle) -> bool {
    let face_data = match brep.faces.get(face) {
        Some(f) => f,
        None => return false,
    };
    let start_he = match brep.loops.get(face_data.outer_loop) {
        Some(l) => l.half_edge,
        None => return false,
    };
    let mut cursor = start_he;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        if brep.half_edges[cursor].origin == vertex {
            return true;
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == start_he {
            break;
        }
    }
    false
}

/// Intersect two 3D line segments. Returns `(t_b, point)` if they cross.
fn line_line_nearest(
    a0: popkern_core::geometry::Point3,
    a1: popkern_core::geometry::Point3,
    b0: popkern_core::geometry::Point3,
    b1: popkern_core::geometry::Point3,
) -> Option<(f64, popkern_core::geometry::Point3)> {
    use popkern_core::geometry::{Point3, Vec3};
    let da = Vec3::new(a1.x - a0.x, a1.y - a0.y, a1.z - a0.z);
    let db = Vec3::new(b1.x - b0.x, b1.y - b0.y, b1.z - b0.z);
    let d0 = Vec3::new(a0.x - b0.x, a0.y - b0.y, a0.z - b0.z);
    let denom = da.dot(&da) * db.dot(&db) - da.dot(&db).powi(2);
    if denom.abs() < 1e-20 {
        return None;
    }
    let t_a = (da.dot(&db) * d0.dot(&db) - db.dot(&db) * d0.dot(&da)) / denom;
    let t_b = (da.dot(&da) * d0.dot(&db) - da.dot(&db) * d0.dot(&da)) / denom;
    if !(-0.1..=1.1).contains(&t_a) || !(-0.1..=1.1).contains(&t_b) {
        return None;
    }
    let pa = Point3::new(a0.x + t_a * da.x, a0.y + t_a * da.y, a0.z + t_a * da.z);
    let pb = Point3::new(b0.x + t_b * db.x, b0.y + t_b * db.y, b0.z + t_b * db.z);
    if pa.distance_to(&pb) > 0.01 {
        return None;
    }
    Some((t_b, pb))
}
