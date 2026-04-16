//! Internal helpers for [`crate::identify`].
//!
//! Holds the boundary-walking, dihedral-angle, and support-pair
//! helpers used by `identify_blend`. Extracted from `identify.rs`
//! so each source file stays under the 500-line structural cap.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Data collected about one boundary half-edge of a candidate face.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BoundaryEntry {
    /// The half-edge on `face`'s outer loop.
    #[allow(dead_code)]
    pub half_edge: HalfEdgeHandle,
    /// The edge beneath `half_edge`.
    #[allow(dead_code)]
    pub edge: EdgeHandle,
    /// The face on the other side of `half_edge` (via the twin).
    pub neighbour: FaceHandle,
    /// The 3D midpoint of the edge — used for sampling surface normals.
    #[allow(dead_code)]
    pub midpoint: Point3,
    /// The dihedral angle between `face`'s normal and `neighbour`'s
    /// normal at `midpoint`, in `[0, π]`. Near 0 = G1 tangent; near
    /// π/2 = sharp 90°; near π = faces pointing anti-parallel.
    pub angle: f64,
}

/// Walk the outer loop of `face` and return one [`BoundaryEntry`]
/// per half-edge. Returns `None` when `face` has no surface, no
/// outer loop, or a broken half-edge chain.
pub(crate) fn boundary_tangency(brep: &BRep, face: FaceHandle) -> Option<Vec<BoundaryEntry>> {
    let face_data = brep.faces.get(face)?;
    let surface_on_face_a = face_data.surface.as_ref()?;
    let surface_a = &surface_on_face_a.geometry;
    let orient_a = surface_on_face_a.orientation;

    let loop_handle = face_data.outer_loop;
    let loop_data = brep.loops.get(loop_handle)?;
    let start = loop_data.half_edge;

    let mut entries = Vec::new();
    let mut he = start;
    let mut safety = 0;
    loop {
        safety += 1;
        if safety > 10_000 {
            // Defensive: broken loop, give up rather than infinite-loop.
            return None;
        }

        let he_data = brep.half_edges.get(he)?;
        let twin_data = brep.half_edges.get(he_data.twin)?;
        let neighbour_loop = brep.loops.get(twin_data.loop_)?;
        let neighbour = neighbour_loop.face;

        let edge = he_data.edge;
        let edge_data = brep.edges.get(edge)?;

        let midpoint = edge_midpoint(brep, edge_data);

        // Skip self-edges (seam): they can't participate in a
        // blend boundary because there's only one face on both sides.
        if neighbour == face {
            he = he_data.next;
            if he == start {
                break;
            }
            continue;
        }

        let Some(neighbour_face_data) = brep.faces.get(neighbour) else {
            he = he_data.next;
            if he == start {
                break;
            }
            continue;
        };
        let Some(surface_on_face_b) = neighbour_face_data.surface.as_ref() else {
            he = he_data.next;
            if he == start {
                break;
            }
            continue;
        };
        let surface_b = &surface_on_face_b.geometry;
        let orient_b = surface_on_face_b.orientation;

        let angle =
            dihedral_angle(surface_a, orient_a, surface_b, orient_b, &midpoint).unwrap_or(f64::NAN);

        if angle.is_finite() {
            entries.push(BoundaryEntry {
                half_edge: he,
                edge,
                neighbour,
                midpoint,
                angle,
            });
        }

        he = he_data.next;
        if he == start {
            break;
        }
    }

    Some(entries)
}

/// Compute the dihedral angle between two faces sharing a 3D point.
///
/// Projects `point` onto each surface to obtain `(u, v)` on both,
/// evaluates the surface normals, and returns the angle between
/// the oriented normals in `[0, π]`.
pub(crate) fn dihedral_angle(
    surface_a: &SurfaceGeometry,
    orient_a: bool,
    surface_b: &SurfaceGeometry,
    orient_b: bool,
    point: &Point3,
) -> Option<f64> {
    let (ua, va) = project_point_to_uv(surface_a, point)?;
    let (ub, vb) = project_point_to_uv(surface_b, point)?;

    let n_a_raw = surface_a.normal_at(ua, va);
    let n_b_raw = surface_b.normal_at(ub, vb);

    // Guard against zero-length normals (parameterisation singularity).
    if n_a_raw.length_squared() < 1e-30 || n_b_raw.length_squared() < 1e-30 {
        return None;
    }

    let n_a = if orient_a { n_a_raw } else { -n_a_raw };
    let n_b = if orient_b { n_b_raw } else { -n_b_raw };

    // Numerically robust angle between unit-ish vectors.
    let cos_theta = (n_a.dot(&n_b) / (n_a.length() * n_b.length())).clamp(-1.0, 1.0);
    Some(cos_theta.acos())
}

/// 3D midpoint of an edge (uses the attached curve if present,
/// otherwise the average of its two endpoint vertex positions).
pub(crate) fn edge_midpoint(brep: &BRep, edge: &popkern_core::topology::Edge) -> Point3 {
    use popkern_core::geometry::curve::CurveEval;

    if let Some(curve_on_edge) = &edge.curve {
        let t_mid = (curve_on_edge.t_start + curve_on_edge.t_end) * 0.5;
        return curve_on_edge.geometry.point_at(t_mid);
    }
    // Fallback: average the vertex positions reached through the
    // half-edge and its twin.
    if let Some(he) = brep.half_edges.get(edge.half_edge)
        && let Some(twin) = brep.half_edges.get(he.twin)
        && let Some(v0) = brep.vertices.get(he.origin)
        && let Some(v1) = brep.vertices.get(twin.origin)
    {
        let a = v0.position;
        let b = v1.position;
        return Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5);
    }
    Point3::new(0.0, 0.0, 0.0)
}

pub(crate) fn distinct_neighbours(boundary: &[BoundaryEntry]) -> Vec<FaceHandle> {
    let mut seen = Vec::new();
    for entry in boundary {
        if !seen.contains(&entry.neighbour) {
            seen.push(entry.neighbour);
        }
    }
    seen
}

/// Extract the first two **tangent** neighbours as the fillet's
/// support-face pair. Returns `None` if fewer than two distinct
/// tangent neighbours exist.
pub(crate) fn tangent_support_pair(
    boundary: &[BoundaryEntry],
    tangent_tol: f64,
) -> Option<(FaceHandle, FaceHandle)> {
    let mut tangent_neighbours: Vec<FaceHandle> = Vec::new();
    for entry in boundary {
        if entry.angle.abs() < tangent_tol && !tangent_neighbours.contains(&entry.neighbour) {
            tangent_neighbours.push(entry.neighbour);
        }
    }
    match tangent_neighbours.len() {
        0 | 1 => None,
        _ => Some((tangent_neighbours[0], tangent_neighbours[1])),
    }
}

/// For a chamfer candidate, return the first two distinct
/// neighbours as the support pair.
pub(crate) fn chamfer_support_pair(
    _boundary: &[BoundaryEntry],
    distinct: &[FaceHandle],
) -> Option<(FaceHandle, FaceHandle)> {
    if distinct.len() < 2 {
        return None;
    }
    Some((distinct[0], distinct[1]))
}

/// Rough chamfer-distance estimate from a chamfer face. For an
/// equal-distance chamfer on a cube edge, the chamfer face is a
/// rectangle whose short side length equals √2 × chamfer-distance.
/// We return half the shortest boundary-edge length as a
/// first-order estimate. Callers that need the exact distance
/// should use metadata (future work) or compute it from the
/// support-face geometry directly.
pub(crate) fn estimate_chamfer_half_width(brep: &BRep, face: FaceHandle) -> Option<f64> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;

    let mut min_len = f64::INFINITY;
    let mut he = start;
    let mut safety = 0;
    loop {
        safety += 1;
        if safety > 10_000 {
            break;
        }
        let he_data = brep.half_edges.get(he)?;
        let edge = brep.edges.get(he_data.edge)?;
        let len = edge_length(brep, edge);
        if len < min_len && len > 0.0 {
            min_len = len;
        }
        he = he_data.next;
        if he == start {
            break;
        }
    }

    if min_len.is_finite() && min_len > 0.0 {
        Some(min_len * 0.5)
    } else {
        None
    }
}

/// Check whether every distinct neighbour of a candidate face has a
/// [`SurfaceGeometry::Plane`] surface.  Used by the Cylinder-as-fillet
/// path: a cylinder face whose neighbours are all planes is a
/// straight-edge plane-plane fillet (§11.1).
pub(crate) fn all_neighbours_are_planes(brep: &BRep, boundary: &[BoundaryEntry]) -> bool {
    let neighbours = distinct_neighbours(boundary);
    if neighbours.is_empty() {
        return false;
    }
    neighbours.iter().all(|&fh| {
        brep.faces
            .get(fh)
            .and_then(|fd| fd.surface.as_ref())
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Plane { .. }))
    })
}

fn edge_length(brep: &BRep, edge: &popkern_core::topology::Edge) -> f64 {
    let Some(he) = brep.half_edges.get(edge.half_edge) else {
        return 0.0;
    };
    let Some(twin) = brep.half_edges.get(he.twin) else {
        return 0.0;
    };
    let Some(v0) = brep.vertices.get(he.origin) else {
        return 0.0;
    };
    let Some(v1) = brep.vertices.get(twin.origin) else {
        return 0.0;
    };
    v0.position.distance_to(&v1.position)
}
