//! Curve clipping to face boundaries.
//!
//! Clips intersection curves to the polygonal boundaries of BRep faces.
//! Used by the intersection stage to trim curves that extend beyond face
//! boundaries due to the UV-box clipping limitations for non-orthogonal
//! plane parameterizations.
//!
//! Line (degree-1) clipping is in the [`clip_line`] submodule.

pub(crate) mod chart_cell;
mod chart_periods;
mod clip_branch_normalize;
mod clip_crossing_reports;
mod clip_crossings_rescue;
mod clip_dispatch;
mod clip_lift_crossings;
mod clip_line;
mod clip_line_exact;
mod clip_negative_branch;
mod clip_pierce_windows;
mod clip_seam_shift;
mod clip_uv;
mod clip_uv_crossings;
mod clip_windows;
mod clip_windows_clearance;
mod clip_windows_decisive;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

pub(crate) use chart_periods::classification_chart_periods;
pub(crate) use clip_dispatch::{clip_curve_to_faces, clip_curve_to_faces_cached};
pub(crate) use clip_line::{clip_curve_to_single_face, clip_line_to_face};
pub(crate) use clip_pierce_windows::PierceInventoryCache;
#[cfg(test)]
use clip_seam_shift::{is_near_seam, shift_curve_u};
#[cfg(test)]
use clip_uv::uv_range_from_crossings;
pub(crate) use clip_uv::{UvCrossing, clip_curve_to_crossings};
pub(crate) use clip_uv_crossings::find_uv_crossings;

/// Stage 2b.58: extract the underlying `CurveOnSurface` (NURBS
/// triple) from an `EdgeTriple` produced by
/// `compute_edge_geometry_unanchored`. Returns `None` if the
/// triple's curves aren't NURBS-backed — shouldn't happen for the
/// marching path the unanchored API uses, but we guard against
/// future analytic-curve paths defensively.
///
/// Bridges the EdgeTriple-shape returned by the new factory with
/// the `CurveOnSurface`-shape that `clip_curve_to_faces` and the
/// rest of `popkern-bool`'s pipeline still consume.
pub(crate) fn edge_triple_to_curve_on_surface(
    triple: &popkern_core::geometry::edge_triple::EdgeTriple,
) -> Option<CurveOnSurface> {
    use popkern_core::geometry::pcurve::PcurveGeometry;
    let CurveGeometry::Nurbs(curve_3d) = &triple.curve_3d else {
        return None;
    };
    let PcurveGeometry::Nurbs(curve_uv1) = &triple.pcurve_a.geometry else {
        return None;
    };
    let PcurveGeometry::Nurbs(curve_uv2) = &triple.pcurve_b.geometry else {
        return None;
    };
    let is_closed = matches!(
        triple.topology,
        popkern_core::geometry::edge_triple::EdgeTopology::Closed { .. }
    );
    Some(CurveOnSurface {
        curve_3d: curve_3d.clone(),
        curve_uv1: curve_uv1.clone(),
        curve_uv2: curve_uv2.clone(),
        is_closed,
    })
}

/// Collect vertex positions from the outer loop of a face by walking half-edges.
///
/// Returns an empty `Vec` if the face's loop has no valid half-edges.
/// Uses a bounded iteration limit to prevent infinite loops on broken topology.
pub(crate) fn face_vertex_positions(brep: &BRep, face_handle: FaceHandle) -> Vec<Point3> {
    let face = &brep.faces[face_handle];
    let loop_handle = face.outer_loop;
    let loop_data = match brep.loops.get(loop_handle) {
        Some(l) => l,
        None => return Vec::new(),
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut positions = Vec::new();
    let mut current = start_he;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = &brep.half_edges[current];
        positions.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}

/// Extract the UV boundary of a face from its half-edge pcurves.
///
/// Returns a list of `(half_edge_handle, pcurve)` forming the outer loop
/// boundary in UV space. Returns `None` if any half-edge lacks a pcurve.
pub(crate) fn face_uv_boundary(
    brep: &BRep,
    face: FaceHandle,
) -> Option<Vec<(HalfEdgeHandle, popkern_core::geometry::pcurve::Pcurve)>> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut boundary = Vec::new();
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = brep.half_edges.get(cursor)?;
        let pcurve = he.pcurve.as_ref()?.clone();
        boundary.push((cursor, pcurve));
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    Some(boundary)
}

/// The UV contours of a face's INNER loops (holes), one walked
/// contour per loop, in the same `(half_edge, pcurve)` form as
/// [`face_uv_boundary`]. A loop whose half-edge chain is broken or
/// lacks a pcurve is skipped rather than failing the face — the
/// outer boundary alone is still a usable (if hole-blind) answer,
/// and the callers that need holes treat these contours additively.
/// Region-semantics consumers (containment parity, record-crossing
/// detection, boundary-clearance hugging) must see holes: a face
/// that gained an inner loop in an earlier boolean otherwise
/// classifies its hole as face interior, and the next boolean
/// imprints trims INSIDE the hole (measured: a lid∪lip union's ring
/// face at the flush plane — the following vent subtract carved its
/// slot into the ring's hole region, OCC `InvalidImbricationOfWires`,
/// 8 open mesh edges).
pub(crate) fn face_uv_inner_boundaries(
    brep: &BRep,
    face: FaceHandle,
) -> Vec<Vec<(HalfEdgeHandle, popkern_core::geometry::pcurve::Pcurve)>> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut contours = Vec::new();
    'loops: for &inner in &face_data.inner_loops {
        let Some(loop_data) = brep.loops.get(inner) else {
            continue;
        };
        let start = loop_data.half_edge;
        if !brep.half_edges.contains_key(start) {
            continue;
        }
        let mut contour = Vec::new();
        let mut cursor = start;
        let limit = brep.half_edges.len() + 1;
        for _ in 0..limit {
            let Some(he) = brep.half_edges.get(cursor) else {
                continue 'loops;
            };
            let Some(pcurve) = he.pcurve.as_ref() else {
                continue 'loops;
            };
            contour.push((cursor, pcurve.clone()));
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
        contours.push(contour);
    }
    contours
}

/// [`face_uv_boundary`] with period branches CHAIN-NORMALIZED: a
/// SEAM edge carries one pcurve shared by both of its uses, so one
/// traversal sits a full period from the walk and the raw contour
/// jumps branches (even-odd containment then reads every interior
/// point as OUTSIDE — a perpendicular bore's essential intersection
/// circle was clipped away against its own cylinder wall and the
/// boolean degraded to the no-records containment path). Each pcurve
/// is shifted by whole periods so its start meets the previous end —
/// the same per-use normalization STEP export performs. Scoped to
/// CONTAINMENT tests; the crossing machinery has its own tuned
/// branch handling and double-shifts if fed normalized input.
#[allow(dead_code)] // prior art kept for the chart-cell campaign
pub(crate) fn face_uv_boundary_chained(
    brep: &BRep,
    face: FaceHandle,
) -> Option<Vec<(HalfEdgeHandle, popkern_core::geometry::pcurve::Pcurve)>> {
    let mut boundary = face_uv_boundary(brep, face)?;
    if let Some((pu, pv)) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| classification_chart_periods(&s.geometry))
        && (pu.is_some() || pv.is_some())
        && boundary.len() > 1
    {
        let mut prev_end: Option<popkern_core::geometry::Point3> = None;
        for (_, pc) in boundary.iter_mut() {
            let (t0, t1) = pc.parameter_range();
            let start_pt = pc.point_at(t0);
            if let Some(prev) = prev_end {
                let mut du = 0.0;
                if let Some(p) = pu {
                    du = ((prev.x - start_pt.x) / p).round() * p;
                }
                let mut dv = 0.0;
                if let Some(p) = pv {
                    dv = ((prev.y - start_pt.y) / p).round() * p;
                }
                if du != 0.0 || dv != 0.0 {
                    *pc = pc.translated_uv(du, dv);
                }
            }
            prev_end = Some(pc.point_at(t1));
        }
    }
    Some(boundary)
}

/// True when EVERY boundary half-edge of `face` has its twin in one
/// of the face's own loops — the face's boundary is nothing but its
/// own seams (a revolved circle's torus: tube seam and equator seam),
/// so it has no genuine boundary and covers its whole chart. Every
/// on-surface point is on such a face.
pub(crate) fn face_covers_full_chart(brep: &BRep, face: FaceHandle) -> bool {
    let Some(boundary) = face_uv_boundary(brep, face) else {
        return false;
    };
    !boundary.is_empty()
        && boundary.iter().all(|(he, _)| {
            brep.half_edges
                .get(*he)
                .and_then(|h| brep.half_edges.get(h.twin))
                .and_then(|t| brep.loops.get(t.loop_))
                .is_some_and(|l| l.face == face)
        })
}

/// PK-304 classification gate: is a torus-surface hit at `(u, v)` on
/// `face`'s chart cell? `Some(true/false)` = decisive; `None` = the
/// hit sits within a border sliver of the cell (ambiguous — belongs
/// to the neighbouring patch as much as this one) or the face is not
/// cell-resolvable, and the caller should treat the ray as poisoned.
pub(crate) fn torus_hit_in_face_cell(
    brep: &BRep,
    face: FaceHandle,
    u: f64,
    v: f64,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Option<bool> {
    let Some(cell) = chart_cell::face_chart_cell(brep, face, tolerance) else {
        // A face with no genuine boundary (both seams of a revolved
        // circle's torus, nothing else) covers its whole chart: every
        // surface hit is on it. Answering "unknown" poisoned every ray
        // from inside the tube and a box subtraction dropped the cut
        // plane's annulus as Outside (PK-418).
        return face_covers_full_chart(brep, face).then_some(true);
    };
    let (fu, fv) = cell.fold(u, v);
    let border = 1.0e-7;
    let near_border = (fu - cell.u_range.0).abs() < border
        || (fu - cell.u_range.1).abs() < border
        || (fv - cell.v_range.0).abs() < border
        || (fv - cell.v_range.1).abs() < border;
    if near_border {
        return None;
    }
    Some(cell.contains(u, v))
}

/// True when any edge of the face's outer loop carries CURVED
/// geometry. The vertex polygon of such a face excludes the bulge
/// between the curve's endpoints, so Cyrus–Beck clips against the
/// CHORD and rejects real trims that live in the bulge (measured: a
/// seat ring's y-side on an ellipse-fronted tilted top was dropped
/// and the split stage failed on the chained remainder). Curved
/// boundaries route to the pcurve-exact UV clip.
pub(crate) fn face_has_curved_boundary(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return false;
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return false;
    }
    let mut he = start;
    let cap = brep.half_edges.len() + 1;
    for _ in 0..cap {
        let d = &brep.half_edges[he];
        if let Some(curve) = brep.edges.get(d.edge).and_then(|e| e.curve.as_ref()) {
            let curved = match &curve.geometry {
                popkern_core::geometry::curve::CurveGeometry::Line { .. } => false,
                popkern_core::geometry::curve::CurveGeometry::Nurbs(c) => c.degree > 1,
                _ => true,
            };
            if curved {
                return true;
            }
        }
        he = d.next;
        if he == start {
            break;
        }
    }
    false
}

/// True when the face's outer vertex polygon turns in both directions
/// — a REFLEX outline. Cyrus–Beck clipping is CONVEX-only: on a
/// reflex polygon its centroid-oriented inward normals disagree and
/// the clip window comes out arbitrary (measured: a seat tool near a
/// reflex deck's crotch kept 300 mm phantom trim lines spanning far
/// outside both bodies; interior trims taking the planar-split
/// fallback fabricated crossings at outline edges). Reflex faces must
/// route to UV clipping.
pub(crate) fn face_polygon_is_reflex(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> bool {
    let pts = brep.face_outer_loop_positions(face);
    if pts.len() < 4 {
        return false;
    }
    // Face normal from Newell's method (robust for any planar
    // polygon), then check turn-direction consistency.
    let n = pts.len();
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    for i in 0..n {
        let p = &pts[i];
        let q = &pts[(i + 1) % n];
        nx += (p.y - q.y) * (p.z + q.z);
        ny += (p.z - q.z) * (p.x + q.x);
        nz += (p.x - q.x) * (p.y + q.y);
    }
    let mut sign = 0.0_f64;
    for i in 0..n {
        let a = &pts[i];
        let b = &pts[(i + 1) % n];
        let c = &pts[(i + 2) % n];
        let e1 = *b - *a;
        let e2 = *c - *b;
        let cross = e1.cross(&e2);
        let s = cross.x * nx + cross.y * ny + cross.z * nz;
        if s.abs() > 1e-9 {
            if sign == 0.0 {
                sign = s.signum();
            } else if s.signum() != sign {
                return true;
            }
        }
    }
    false
}

/// A trim residue emitted by the clip: owner side (`true` = face A),
/// the single-face splitting curve, and per-end flags marking which
/// ends abut a kept sibling's cut bound (PK-472).
pub(crate) type TrimResidue = (bool, CurveOnSurface, (bool, bool));

#[cfg(test)]
#[path = "clip_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "clip_uv_tests.rs"]
mod uv_tests;
