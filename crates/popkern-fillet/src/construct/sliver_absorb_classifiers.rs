//! Classifier helpers for the sliver absorption pass.
//!
//! Extracted from `sliver_absorb.rs` to keep that file under the
//! 500-line structural cap. These routines scan the BRep for a
//! single sliver of each class and return the half-edge that
//! [`crate::construct::sliver_absorb::absorb_slivers`] should pass
//! to `BRep::kef`. They never mutate the BRep — mutation is the
//! driver loop's responsibility.
//!
//! Visibility: items here are `pub(super)` so the parent
//! `sliver_absorb` module can call them and re-export via `use`
//! for its `#[path]`-mounted test module.
//!
//! See the module-level docs of `sliver_absorb.rs` for the
//! architectural invariants (restart-after-mutation, directional
//! kef, iteration cap).
//!
//! No function bodies were modified during the extraction — these
//! are byte-for-byte the originals from `sliver_absorb.rs`.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Class D: ≤4-vert Plane face whose outer loop contains a
/// half-edge whose twin face is a Sphere or Gregory cap. Returns
/// `(plane_face, cap_side_half_edge)` so `kef` preserves the
/// sphere labelling.
pub(super) fn find_corner_fragment_target(
    brep: &BRep,
    retained: &std::collections::HashSet<FaceHandle>,
    legacy_caps: &std::collections::BTreeSet<FaceHandle>,
) -> Option<(FaceHandle, HalfEdgeHandle)> {
    for fh in brep.faces.keys() {
        if retained.contains(&fh) {
            continue;
        }
        let face = &brep.faces[fh];
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
            continue;
        }
        let positions = brep.face_outer_loop_positions(fh);
        if positions.len() > 4 {
            continue;
        }
        let start_he = brep.loops[face.outer_loop].half_edge;
        let mut cursor = start_he;
        let walk_cap = brep.half_edges.len() + 1;
        for _ in 0..walk_cap {
            let he = &brep.half_edges[cursor];
            let twin_h = he.twin;
            let Some(twin) = brep.half_edges.get(twin_h) else {
                cursor = he.next;
                continue;
            };
            let Some(twin_loop) = brep.loops.get(twin.loop_) else {
                cursor = he.next;
                continue;
            };
            let twin_face = twin_loop.face;
            // Guarded like the loop lookup above: mid-pipeline
            // topology can carry a stale loop→face reference after a
            // blend surgery (measured: a chamfer on a slotted plate
            // panicked here on an invalid SlotMap key).
            // A declared legacy vertex-blend cap absorbs its corner
            // fragments regardless of surface kind: conic/variable
            // profile caps assign a plain fitted NURBS (PK-222), and
            // gating on the surface enum starved the whole downstream
            // welded-corner reconstruction for every non-sphere,
            // non-Gregory cap.
            if twin_face != fh
                && let Some(twin_face_data) = brep.faces.get(twin_face)
                && let Some(twin_sof) = twin_face_data.surface.as_ref()
                && (matches!(
                    twin_sof.geometry,
                    SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Gregory(_)
                ) || legacy_caps.contains(&twin_face))
            {
                // A corner FRAGMENT is commensurate with the cap it
                // merges into; without a size gate this class ate
                // entire trimmed side walls (≤4 vertices after
                // trimming, adjacent to a sphere cap) and cascaded
                // 16 absorptions on a boss-plate chain fillet. Gate
                // on polygonal area relative to the cap's own.
                let poly_area = |handle: FaceHandle| -> f64 {
                    let pts = brep.face_outer_loop_positions(handle);
                    if pts.len() < 3 {
                        return 0.0;
                    }
                    let p0 = pts[0];
                    let mut total = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
                    for i in 1..pts.len() - 1 {
                        let v1 = pts[i] - p0;
                        let v2 = pts[i + 1] - p0;
                        total = total + v1.cross(&v2);
                    }
                    total.length() * 0.5
                };
                let frag_area = poly_area(fh);
                let cap_area = poly_area(twin_face);
                if frag_area <= cap_area * 4.0 {
                    return Some((fh, twin_h));
                }
            }
            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }
    }
    None
}

/// Class A: face whose outer-loop polygonal area is below
/// `tolerance.linear()²`. Returns `(sliver_face, neighbour_side_he, area)`.
pub(super) fn find_zero_area_target(
    brep: &BRep,
    tolerance: &Tolerance,
    retained: &std::collections::HashSet<FaceHandle>,
) -> Option<(FaceHandle, HalfEdgeHandle, f64)> {
    let tol_sq = tolerance.linear() * tolerance.linear();
    for fh in brep.faces.keys() {
        if retained.contains(&fh) {
            continue;
        }
        let positions = brep.face_outer_loop_positions(fh);
        if positions.len() < 3 {
            continue;
        }
        let p0 = positions[0];
        let mut total = Vec3::new(0.0, 0.0, 0.0);
        for i in 1..positions.len() - 1 {
            let v1 = positions[i] - p0;
            let v2 = positions[i + 1] - p0;
            total = total + v1.cross(&v2);
        }
        let area = total.length() * 0.5;
        if area >= tol_sq {
            continue;
        }
        // Zero-area face. Find a half-edge on its outer loop whose
        // twin face differs — return the TWIN half-edge so the
        // neighbour survives the merge.
        let face_data = &brep.faces[fh];
        let start_he = brep.loops[face_data.outer_loop].half_edge;
        let mut cursor = start_he;
        let walk_cap = brep.half_edges.len() + 1;
        for _ in 0..walk_cap {
            let he = &brep.half_edges[cursor];
            let twin_h = he.twin;
            let Some(twin) = brep.half_edges.get(twin_h) else {
                cursor = he.next;
                continue;
            };
            let Some(twin_loop) = brep.loops.get(twin.loop_) else {
                cursor = he.next;
                continue;
            };
            let twin_face = twin_loop.face;
            if twin_face != fh {
                return Some((fh, twin_h, area));
            }
            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }
    }
    None
}

/// Class C: two adjacent faces share an edge, and their underlying
/// surfaces are structurally equal (same type, same params within
/// tolerance). Return the half-edge of the face with MORE outer-
/// loop vertices so the larger face survives the merge.
pub(super) fn find_redundant_edge_target(
    brep: &BRep,
    tolerance: &Tolerance,
    retained: &std::collections::HashSet<FaceHandle>,
) -> Option<(FaceHandle, HalfEdgeHandle)> {
    for eh in brep.edges.keys() {
        let he_h = brep.edges[eh].half_edge;
        let Some(he) = brep.half_edges.get(he_h) else {
            continue;
        };
        let twin_h = he.twin;
        let Some(twin) = brep.half_edges.get(twin_h) else {
            continue; // boundary edge
        };
        let Some(loop_a) = brep.loops.get(he.loop_) else {
            continue;
        };
        let Some(loop_b) = brep.loops.get(twin.loop_) else {
            continue;
        };
        let face_a = loop_a.face;
        let face_b = loop_b.face;
        if face_a == face_b {
            continue; // seam edge
        }
        if retained.contains(&face_a) || retained.contains(&face_b) {
            continue;
        }
        let Some(sa) = brep
            .faces
            .get(face_a)
            .and_then(|face| face.surface.as_ref())
        else {
            continue;
        };
        let Some(sb) = brep
            .faces
            .get(face_b)
            .and_then(|face| face.surface.as_ref())
        else {
            continue;
        };
        if !surfaces_equal(&sa.geometry, &sb.geometry, tolerance) {
            continue;
        }
        let verts_a = brep.face_outer_loop_positions(face_a).len();
        let verts_b = brep.face_outer_loop_positions(face_b).len();
        // kef kills the TWIN face and keeps the half-edge's face.
        // Return the half-edge on the larger face so the smaller
        // one dies in the merge.
        if verts_a >= verts_b {
            return Some((face_b, he_h));
        } else {
            return Some((face_a, twin_h));
        }
    }
    None
}

/// Class B: face whose outer loop contains an edge shorter than
/// `tolerance.linear()`. Returns `(sliver_face, short_edge_he, min_length)`.
pub(super) fn find_short_edge_target(
    brep: &BRep,
    tolerance: &Tolerance,
    retained: &std::collections::HashSet<FaceHandle>,
) -> Option<(FaceHandle, HalfEdgeHandle, f64)> {
    let tol = tolerance.linear();
    for fh in brep.faces.keys() {
        if retained.contains(&fh) {
            continue;
        }
        let Some(face_data) = brep.faces.get(fh) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start_he = loop_data.half_edge;
        let mut cursor = start_he;
        let walk_cap = brep.half_edges.len() + 1;
        for _ in 0..walk_cap {
            let Some(he) = brep.half_edges.get(cursor) else {
                break;
            };
            let Some(next_he) = brep.half_edges.get(he.next) else {
                break;
            };
            let Some(start_vertex) = brep.vertices.get(he.origin) else {
                break;
            };
            let Some(end_vertex) = brep.vertices.get(next_he.origin) else {
                break;
            };
            let p_start = start_vertex.position;
            let p_end = end_vertex.position;
            let length = (p_end - p_start).length();
            if length < tol {
                // Short edge. Return `cursor` — the half-edge
                // whose face is the current face. kef will kill
                // the neighbour across this edge, which may or
                // may not be the right choice, so the caller
                // relies on `retained_faces` to bail out if the
                // kef fails.
                return Some((fh, cursor, length));
            }
            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }
    }
    None
}

/// Structural equality check for two `SurfaceGeometry` values
/// within tolerance. Used by Class C (redundant-edge) to decide
/// whether two adjacent faces can be merged without changing the
/// geometric meaning of the result.
///
/// Covers the analytic primitives (Plane, Cylinder, Sphere, Cone,
/// Torus). `Nurbs` comparisons return `false` — control-polygon
/// equality is expensive and rarely meaningful for fillet output
/// anyway.
pub(super) fn surfaces_equal(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> bool {
    use popkern_core::geometry::Point3;
    let tol = tolerance.linear();
    let tol_ang = tolerance.angular();
    let point_close = |p: &Point3, q: &Point3| -> bool { p.distance_to(q) < tol };
    let dir_parallel = |u: &Vec3, v: &Vec3| -> bool {
        let nu = u.normalized().unwrap_or(Vec3::new(0.0, 0.0, 0.0));
        let nv = v.normalized().unwrap_or(Vec3::new(0.0, 0.0, 0.0));
        (nu.dot(&nv).abs() - 1.0).abs() < tol_ang
    };
    match (a, b) {
        (
            SurfaceGeometry::Plane {
                origin: oa,
                normal: na,
                ..
            },
            SurfaceGeometry::Plane {
                origin: ob,
                normal: nb,
                ..
            },
        ) => {
            // Planes are equal if normals are parallel (same or
            // opposite direction) and the origins lie in the same
            // plane (offset along the normal is zero).
            if !dir_parallel(na, nb) {
                return false;
            }
            let delta = *ob - *oa;
            delta.dot(na).abs() < tol
        }
        (
            SurfaceGeometry::Cylinder {
                origin: oa,
                axis: aa,
                radius: ra,
            },
            SurfaceGeometry::Cylinder {
                origin: ob,
                axis: ab,
                radius: rb,
            },
        ) => {
            if !dir_parallel(aa, ab) {
                return false;
            }
            if (ra - rb).abs() >= tol {
                return false;
            }
            // Origins must lie on the same axis line.
            let delta = *ob - *oa;
            let along = delta.dot(&aa.normalized().unwrap_or(Vec3::new(0.0, 0.0, 0.0)));
            let perp = delta - aa.normalized().unwrap_or(Vec3::new(0.0, 0.0, 0.0)) * along;
            perp.length() < tol
        }
        (
            SurfaceGeometry::Sphere {
                center: ca,
                radius: ra,
            },
            SurfaceGeometry::Sphere {
                center: cb,
                radius: rb,
            },
        ) => point_close(ca, cb) && (ra - rb).abs() < tol,
        (
            SurfaceGeometry::Cone {
                apex: aa,
                axis: ax_a,
                half_angle: ha,
            },
            SurfaceGeometry::Cone {
                apex: ab,
                axis: ax_b,
                half_angle: hb,
            },
        ) => point_close(aa, ab) && dir_parallel(ax_a, ax_b) && (ha - hb).abs() < tol_ang,
        (
            SurfaceGeometry::Torus {
                center: ca,
                axis: aa,
                major_radius: ma,
                minor_radius: na,
            },
            SurfaceGeometry::Torus {
                center: cb,
                axis: ab,
                major_radius: mb,
                minor_radius: nb,
            },
        ) => {
            point_close(ca, cb)
                && dir_parallel(aa, ab)
                && (ma - mb).abs() < tol
                && (na - nb).abs() < tol
        }
        // Mismatched surface types are never equal. NURBS vs NURBS
        // would need control-polygon comparison which is rarely
        // useful here — leave `false`.
        // NURBS: structural identity (cloned assignments from the
        // operation's fillet-surface list — the vertical band strips
        // carried eg NURBS surfaces that this equality couldn't see,
        // so surface_fix's sphere fallback ate them).
        (SurfaceGeometry::Nurbs(na), SurfaceGeometry::Nurbs(nb)) => {
            na.degree_u == nb.degree_u
                && na.degree_v == nb.degree_v
                && na.num_u == nb.num_u
                && na.num_v == nb.num_v
                && na.knots_u == nb.knots_u
                && na.knots_v == nb.knots_v
                && na.control_points == nb.control_points
        }
        _ => false,
    }
}
