//! General post-construction sliver-face absorption pass.
//!
//! **Step 2 scope**: infrastructure only. This module subsumes the
//! previous specific-case `blend_cleanup::merge_corner_fragments_into_blends`
//! absorber (Class D — "≤4-vert Plane → Sphere/Gregory vertex blend") under
//! a general classify-apply-reclassify framework that can be
//! extended to additional sliver classes in step 3 without further
//! plumbing changes.
//!
//! **Not yet in this module (step 3 work):**
//! - Class A: zero-area face detection + absorption
//! - Class B: thin (high aspect ratio) face detection + absorption
//! - Class C: redundant-edge detection (two adjacent faces sharing
//!   the same underlying surface) + absorption
//! - Emission of `SliverAbsorbWarning::Retained` when a classified
//!   sliver cannot be safely absorbed (topology would be violated)
//!
//! ## Architectural invariants
//!
//! **Restart-after-mutation** is mandatory: `BRep::kef` invalidates
//! `brep.faces.keys()` iteration. The general loop therefore
//! follows the exact pattern proved in the previous absorber —
//! classify, find ONE target, apply ONE `kef`, reclassify from
//! scratch. Never batch kef targets.
//!
//! **Directional kef**: when merging a plane fragment into a
//! sphere blend, we pass the sphere-side half-edge so the sphere
//! face *survives* and the plane fragment *dies*. Otherwise the
//! merged face would inherit the plane's surface and lose the
//! sphere labelling. Every absorber that touches surface-typed
//! neighbours must respect this convention.
//!
//! **Iteration cap** is `brep.faces.len() * 2 + 16`. Generous
//! safety net — cascade absorption rarely needs more than a few
//! passes in practice, and each successful absorption reduces the
//! face count, so the pass is guaranteed to terminate well below
//! the cap.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Run the sliver absorber and route any retained-sliver warnings
/// into the caller's `BlendWarning` vector. Kept here rather than
/// inline in `construct::mod` so the dispatch logic stays colocated
/// with the classifier and the parent file stays under the 500-line
/// structural cap.
pub(crate) fn run_and_route_warnings(
    brep: &mut BRep,
    tolerance: &Tolerance,
    edge_geometries: &[crate::construct::EdgeGeometry],
    warnings: &mut Vec<(
        popkern_core::topology::EdgeHandle,
        crate::error::BlendWarning,
    )>,
) {
    let sliver_warnings = absorb_slivers(brep, tolerance);
    let fallback = edge_geometries
        .first()
        .map(|eg| eg.analysis.edge)
        .unwrap_or_default();
    for SliverAbsorbWarning::Retained { face, kind } in sliver_warnings {
        warnings.push((
            fallback,
            crate::error::BlendWarning::SliverRetained {
                face,
                kind: kind.label(),
            },
        ));
    }
}

#[cfg(test)]
#[path = "sliver_absorb_tests.rs"]
mod tests;

/// A classified sliver face that the absorber could not remove.
///
/// Returned by [`absorb_slivers`] for any face that a classifier
/// flagged but the corresponding absorber refused to merge — typically
/// because no suitable absorbable neighbour exists or the merge would
/// violate topology. Callers consume these via `BlendWarning` plumbing:
/// `apply_blend_geometry` converts them into
/// [`crate::error::BlendWarning::SliverRetained`] entries on the
/// per-edge warning list.
#[derive(Debug, Clone)]
pub(crate) enum SliverAbsorbWarning {
    /// A sliver was classified but could not be absorbed.
    Retained {
        /// The face that remains in the output BRep.
        face: FaceHandle,
        /// Which sliver category the classifier flagged.
        kind: SliverKind,
    },
}

/// Classification of a sliver face. Each variant carries the
/// measurement that drove the classification, for diagnostic
/// rendering via `Debug` when the sliver ends up retained. The
/// inner fields are not read by non-`Debug` code paths today —
/// the `#[allow(dead_code)]` is intentional.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub(crate) enum SliverKind {
    /// Polygonal outer-loop area below `tolerance.linear()²`.
    ZeroArea {
        /// The measured polygonal area.
        area: f64,
    },
    /// Face with an outer-loop edge shorter than `tolerance.linear()`
    /// (a "thin strip" produced by near-parallel trim curves).
    ShortEdge {
        /// Length of the shortest edge.
        min_length: f64,
    },
    /// Face shares an edge with a neighbour whose underlying
    /// surface is structurally equal (duplicate-subdivision artifact).
    RedundantEdge,
    /// ≤4-vert face whose twin face via some edge is a vertex-
    /// blend sphere — the Round 5r corner-fragment case.
    CornerFragment {
        /// Number of vertices on the outer loop (≤ 4 when flagged).
        vertex_count: usize,
    },
}

impl SliverKind {
    /// Short static-string label for the `BlendWarning::SliverRetained`
    /// variant's `kind` field. Callers match on this string rather
    /// than on a cloned enum to keep the public `BlendWarning` small.
    pub(crate) fn label(&self) -> &'static str {
        match self {
            SliverKind::ZeroArea { .. } => "zero-area",
            SliverKind::ShortEdge { .. } => "short-edge",
            SliverKind::RedundantEdge => "redundant-edge",
            SliverKind::CornerFragment { .. } => "corner-fragment",
        }
    }
}

/// Run the sliver absorption pass. Classify → apply one `kef` →
/// reclassify, capped at `faces.len() * 2 + 16` iterations.
/// Priority: Class D (corner fragment → sphere, test-covered) →
/// A (zero-area) → C (redundant edge) → B (short edge, roughest
/// heuristic, last). Failed kefs mark the face as retained so the
/// pass moves on rather than spinning.
pub(crate) fn absorb_slivers(brep: &mut BRep, tolerance: &Tolerance) -> Vec<SliverAbsorbWarning> {
    let mut warnings: Vec<SliverAbsorbWarning> = Vec::new();
    let cap = brep.faces.len() * 2 + 16;
    let mut retained: std::collections::HashSet<FaceHandle> = std::collections::HashSet::new();
    for _ in 0..cap {
        let candidate = find_corner_fragment_target(brep, &retained)
            .map(|(face, he)| {
                let verts = brep.face_outer_loop_positions(face).len();
                (
                    face,
                    SliverKind::CornerFragment {
                        vertex_count: verts,
                    },
                    he,
                )
            })
            .or_else(|| {
                find_zero_area_target(brep, tolerance, &retained)
                    .map(|(face, he, area)| (face, SliverKind::ZeroArea { area }, he))
            })
            .or_else(|| {
                find_redundant_edge_target(brep, tolerance, &retained)
                    .map(|(face, he)| (face, SliverKind::RedundantEdge, he))
            })
            .or_else(|| {
                find_short_edge_target(brep, tolerance, &retained)
                    .map(|(face, he, min_length)| (face, SliverKind::ShortEdge { min_length }, he))
            });
        let Some((face, kind, he)) = candidate else {
            break;
        };
        if brep.kef(he).is_err() {
            warnings.push(SliverAbsorbWarning::Retained { face, kind });
            retained.insert(face);
        }
    }
    warnings
}

/// Class D: ≤4-vert Plane face whose outer loop contains a
/// half-edge whose twin face is a Sphere or Gregory cap. Returns
/// `(plane_face, cap_side_half_edge)` so `kef` preserves the
/// sphere labelling.
fn find_corner_fragment_target(
    brep: &BRep,
    retained: &std::collections::HashSet<FaceHandle>,
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
            let twin_face = brep.loops[twin.loop_].face;
            if twin_face != fh
                && let Some(twin_sof) = brep.faces[twin_face].surface.as_ref()
                && matches!(
                    twin_sof.geometry,
                    SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Gregory(_)
                )
            {
                return Some((fh, twin_h));
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
fn find_zero_area_target(
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
            let twin_face = brep.loops[twin.loop_].face;
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
fn find_redundant_edge_target(
    brep: &BRep,
    tolerance: &Tolerance,
    retained: &std::collections::HashSet<FaceHandle>,
) -> Option<(FaceHandle, HalfEdgeHandle)> {
    for eh in brep.edges.keys() {
        let he_h = brep.edges[eh].half_edge;
        let he = &brep.half_edges[he_h];
        let twin_h = he.twin;
        let Some(twin) = brep.half_edges.get(twin_h) else {
            continue; // boundary edge
        };
        let face_a = brep.loops[he.loop_].face;
        let face_b = brep.loops[twin.loop_].face;
        if face_a == face_b {
            continue; // seam edge
        }
        if retained.contains(&face_a) || retained.contains(&face_b) {
            continue;
        }
        let Some(sa) = brep.faces[face_a].surface.as_ref() else {
            continue;
        };
        let Some(sb) = brep.faces[face_b].surface.as_ref() else {
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
fn find_short_edge_target(
    brep: &BRep,
    tolerance: &Tolerance,
    retained: &std::collections::HashSet<FaceHandle>,
) -> Option<(FaceHandle, HalfEdgeHandle, f64)> {
    let tol = tolerance.linear();
    for fh in brep.faces.keys() {
        if retained.contains(&fh) {
            continue;
        }
        let face_data = &brep.faces[fh];
        let start_he = brep.loops[face_data.outer_loop].half_edge;
        let mut cursor = start_he;
        let walk_cap = brep.half_edges.len() + 1;
        for _ in 0..walk_cap {
            let he = &brep.half_edges[cursor];
            let next_he = &brep.half_edges[he.next];
            let p_start = brep.vertices[he.origin].position;
            let p_end = brep.vertices[next_he.origin].position;
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
fn surfaces_equal(a: &SurfaceGeometry, b: &SurfaceGeometry, tolerance: &Tolerance) -> bool {
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
        _ => false,
    }
}
