//! Detector implementations for `popkern_core::diagnostics`.
//!
//! Split out via `#[path]` to keep the public-API file under the
//! 500-line ceiling. Each detector is a free function returning a
//! `Vec<DefectKind>`; `Diagnostics::run` orchestrates them.

use std::collections::BTreeMap;

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, FaceHandle, VertexHandle};

#[path = "diagnostics_helpers.rs"]
mod helpers;
use helpers::{
    aabb_overlap, check_nurbs_curve, check_nurbs_surface, edge_length, edge_midpoint, face_area,
    face_bbox, face_centroid, face_edge_aspect_ratio, min_median_max, surface_kind_label,
    surfaces_match,
};

use super::{
    CoincidentPair, Defects, DegenerateEdge, FaceOverlap, GeometryStats, NakedEdge, NurbsDefect,
    NurbsKind, NurbsScope, PcurveEscape, SliverFace, TJunction, TopologyStats,
};

// ---------------------------------------------------------------------------
// Aggregate
// ---------------------------------------------------------------------------

pub(super) fn topology_stats(brep: &BRep) -> TopologyStats {
    let mut surface_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut faces_with_surface = 0;
    for (_, face) in &brep.faces {
        if let Some(sof) = face.surface.as_ref() {
            faces_with_surface += 1;
            let key = surface_kind_label(&sof.geometry);
            *surface_kinds.entry(key).or_insert(0) += 1;
        }
    }
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    TopologyStats {
        faces: brep.faces.len(),
        faces_with_surface,
        edges: brep.edges.len(),
        vertices: brep.vertices.len(),
        half_edges: brep.half_edges.len(),
        euler_chi: v - e + f,
        surface_kinds,
    }
}

pub(super) fn geometry_stats(brep: &BRep) -> GeometryStats {
    let edge_lengths: Vec<f64> = brep.edges.keys().map(|eh| edge_length(brep, eh)).collect();
    let face_areas: Vec<f64> = brep.faces.keys().map(|fh| face_area(brep, fh)).collect();
    let (e_min, e_med, e_max) = min_median_max(&edge_lengths);
    let (f_min, f_med, f_max) = min_median_max(&face_areas);
    let (bbox_min, bbox_max) = vertex_bbox(brep);
    let dx = bbox_max.x - bbox_min.x;
    let dy = bbox_max.y - bbox_min.y;
    let dz = bbox_max.z - bbox_min.z;
    GeometryStats {
        edge_length_min: e_min,
        edge_length_median: e_med,
        edge_length_max: e_max,
        face_area_min: f_min,
        face_area_median: f_med,
        face_area_max: f_max,
        bbox_min,
        bbox_max,
        bbox_diagonal: (dx * dx + dy * dy + dz * dz).sqrt(),
    }
}

fn vertex_bbox(brep: &BRep) -> (Point3, Point3) {
    let mut min = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any = false;
    for (_, v) in &brep.vertices {
        any = true;
        let p = v.position;
        if p.x < min.x {
            min.x = p.x;
        }
        if p.y < min.y {
            min.y = p.y;
        }
        if p.z < min.z {
            min.z = p.z;
        }
        if p.x > max.x {
            max.x = p.x;
        }
        if p.y > max.y {
            max.y = p.y;
        }
        if p.z > max.z {
            max.z = p.z;
        }
    }
    if !any {
        return (Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0));
    }
    (min, max)
}

pub(super) fn detect_defects(brep: &BRep) -> Defects {
    let (count, max_div) = tess_non_manifold_summary(brep);
    Defects {
        naked_edges: detect_naked_edges(brep),
        degenerate_edges: detect_degenerate_edges(brep),
        sliver_faces: detect_sliver_faces(brep),
        coincident_vertices: detect_coincident_vertices(brep),
        t_junctions: detect_t_junctions(brep),
        nurbs_defects: detect_nurbs_defects(brep),
        pcurve_escapes: detect_pcurve_escapes(brep),
        face_overlaps: detect_face_overlaps(brep),
        tess_non_manifold_count: count,
        tess_non_manifold_max: max_div,
    }
}

fn tess_non_manifold_summary(brep: &BRep) -> (usize, f64) {
    let threshold = 5.0 * brep.tolerance.linear().max(1e-9);
    let reports = brep.find_tessellation_non_manifold_edges(16, threshold);
    let max_div = reports
        .iter()
        .map(|r| r.max_divergence)
        .fold(0.0_f64, f64::max);
    (reports.len(), max_div)
}

// ---------------------------------------------------------------------------
// Round 1: naked, degenerate, sliver, coincident
// ---------------------------------------------------------------------------

fn detect_naked_edges(brep: &BRep) -> Vec<NakedEdge> {
    let mut out = Vec::new();
    for (eh, edge) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let Some(loop_a) = brep.loops.get(he.loop_) else {
            continue;
        };
        let Some(loop_b) = brep.loops.get(twin.loop_) else {
            continue;
        };
        let face_a = loop_a.face;
        let face_b = loop_b.face;
        let face_a_has = brep.faces.get(face_a).is_some_and(|f| f.surface.is_some());
        let face_b_has = brep.faces.get(face_b).is_some_and(|f| f.surface.is_some());
        if face_a_has && face_b_has {
            continue;
        }
        out.push(NakedEdge {
            edge: eh,
            face_a,
            face_a_has_surface: face_a_has,
            face_b,
            face_b_has_surface: face_b_has,
            midpoint: edge_midpoint(brep, eh),
        });
    }
    out
}

fn detect_degenerate_edges(brep: &BRep) -> Vec<DegenerateEdge> {
    let tol = brep.tolerance.linear().max(1e-12);
    let mut out = Vec::new();
    for (eh, _) in &brep.edges {
        let length = edge_length(brep, eh);
        if length < tol {
            out.push(DegenerateEdge {
                edge: eh,
                length,
                midpoint: edge_midpoint(brep, eh),
            });
        }
    }
    out
}

fn detect_sliver_faces(brep: &BRep) -> Vec<SliverFace> {
    let tol = brep.tolerance.linear().max(1e-12);
    let area_floor = tol * tol * 100.0;
    let aspect_ceiling = 1000.0;
    let mut out = Vec::new();
    for (fh, _) in &brep.faces {
        let area = face_area(brep, fh);
        let aspect = face_edge_aspect_ratio(brep, fh);
        if area < area_floor || aspect > aspect_ceiling {
            out.push(SliverFace {
                face: fh,
                area,
                aspect_ratio: aspect,
                centroid: face_centroid(brep, fh),
            });
        }
    }
    out
}

fn detect_coincident_vertices(brep: &BRep) -> Vec<CoincidentPair> {
    let tol = brep.tolerance.linear().max(1e-12);
    let handles: Vec<VertexHandle> = brep.vertices.keys().collect();
    let mut out = Vec::new();
    for i in 0..handles.len() {
        let a = handles[i];
        let pa = brep.vertices[a].position;
        for &b in &handles[i + 1..] {
            let pb = brep.vertices[b].position;
            let d = pa.distance_to(&pb);
            if d < tol {
                out.push(CoincidentPair {
                    a,
                    b,
                    distance: d,
                    midpoint: Point3::new(
                        0.5 * (pa.x + pb.x),
                        0.5 * (pa.y + pb.y),
                        0.5 * (pa.z + pb.z),
                    ),
                });
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Round 2: T-junctions, NURBS irregularities, pcurve UV-escape
// ---------------------------------------------------------------------------

/// Detect vertices that lie on the *interior* of an edge curve they
/// aren't an endpoint of — these break manifold-edge invariants and
/// confuse OCCT-strict slicers. Brute force O(V × E) since BRep
/// fixtures are small (tens of edges, tens of vertices).
fn detect_t_junctions(brep: &BRep) -> Vec<TJunction> {
    let tol = brep.tolerance.linear().max(1e-12);
    let on_edge_threshold = tol * 5.0;
    let mut out = Vec::new();
    for (eh, edge) in &brep.edges {
        let Some(curve_on_edge) = edge.curve.as_ref() else {
            continue;
        };
        let (he_origin, twin_origin) = match (
            brep.half_edges.get(edge.half_edge),
            brep.half_edges
                .get(edge.half_edge)
                .and_then(|h| brep.half_edges.get(h.twin)),
        ) {
            (Some(h), Some(t)) => (h.origin, t.origin),
            _ => continue,
        };
        let edge_len = edge_length(brep, eh);
        let interior_param_margin = (curve_on_edge.t_end - curve_on_edge.t_start) * 0.05;
        for (vh, v) in &brep.vertices {
            if vh == he_origin || vh == twin_origin {
                continue;
            }
            let p = v.position;
            let Some(t) = curve_on_edge.geometry.closest_point(&p) else {
                continue;
            };
            // Must be strictly interior — skip projections that
            // fall on (or near) the endpoints.
            if t < curve_on_edge.t_start + interior_param_margin
                || t > curve_on_edge.t_end - interior_param_margin
            {
                continue;
            }
            let p_on_curve = curve_on_edge.geometry.point_at(t);
            let d = p.distance_to(&p_on_curve);
            if d < on_edge_threshold && d < edge_len * 0.5 {
                out.push(TJunction {
                    vertex: vh,
                    edge: eh,
                    parameter: t,
                    distance: d,
                    vertex_position: p,
                });
            }
        }
    }
    out
}

/// Detect ill-formed NURBS curves and surfaces in the BRep — invalid
/// weights, knot multiplicity exceeding degree+1, mismatched control-
/// point/knot counts, non-finite control points, or non-monotonic knots.
/// These are upstream bugs in the NURBS construction code.
fn detect_nurbs_defects(brep: &BRep) -> Vec<NurbsDefect> {
    let mut out = Vec::new();
    // Edge curves
    for (eh, edge) in &brep.edges {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        if let crate::geometry::curve::CurveGeometry::Nurbs(nurbs) = &curve.geometry {
            for kind in check_nurbs_curve(nurbs) {
                out.push(NurbsDefect {
                    scope: NurbsScope::EdgeCurve(eh),
                    kind,
                });
            }
        }
    }
    // Half-edge pcurves — analytic variants (Linear) are inherently
    // valid (no CPs to check); only the NURBS form goes through the
    // NURBS defect checker.
    for (heh, he) in &brep.half_edges {
        if let Some(pc) = he.pcurve.as_ref()
            && let Some(nurbs) = pc.as_nurbs()
        {
            for kind in check_nurbs_curve(nurbs) {
                out.push(NurbsDefect {
                    scope: NurbsScope::HalfEdgePcurve(heh),
                    kind,
                });
            }
        }
    }
    // Face surfaces — Round 3 extends NURBS validation to surfaces.
    for (fh, face) in &brep.faces {
        if let Some(sof) = face.surface.as_ref()
            && let SurfaceGeometry::Nurbs(n) = &sof.geometry
        {
            for kind in check_nurbs_surface(n) {
                out.push(NurbsDefect {
                    scope: NurbsScope::FaceSurface(fh),
                    kind,
                });
            }
        }
    }
    out
}

/// Detect pairs of faces that look like duplicates: same surface
/// kind with matching parameters AND overlapping vertex AABBs AND
/// no shared edge. O(F²) brute force — fine for tens of faces.
///
/// Adjacent faces (e.g. corner-blend caps at a pyramid apex)
/// trivially have overlapping AABBs because they share an edge —
/// those are *not* duplicates. Real duplicate / un-merged faces
/// from sliver_absorb or boolean ops cover the same 3-D region
/// without sharing edge topology.
fn detect_face_overlaps(brep: &BRep) -> Vec<FaceOverlap> {
    let tol = brep.tolerance.linear().max(1e-12);
    let handles: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut bboxes: Vec<((Point3, Point3), Point3)> = Vec::with_capacity(handles.len());
    let mut face_edges: Vec<std::collections::BTreeSet<crate::topology::EdgeHandle>> =
        Vec::with_capacity(handles.len());
    for &fh in &handles {
        bboxes.push((face_bbox(brep, fh), face_centroid(brep, fh)));
        face_edges.push(face_edge_set(brep, fh));
    }
    let mut out = Vec::new();
    for i in 0..handles.len() {
        for j in (i + 1)..handles.len() {
            let fa = handles[i];
            let fb = handles[j];
            if !aabb_overlap(&bboxes[i].0, &bboxes[j].0, tol * 2.0) {
                continue;
            }
            // Skip pairs sharing any edge — they're adjacent, not
            // duplicates.
            if face_edges[i].intersection(&face_edges[j]).next().is_some() {
                continue;
            }
            let surf_a = brep.faces.get(fa).and_then(|f| f.surface.as_ref());
            let surf_b = brep.faces.get(fb).and_then(|f| f.surface.as_ref());
            let (Some(sa), Some(sb)) = (surf_a, surf_b) else {
                continue;
            };
            if !surfaces_match(&sa.geometry, &sb.geometry, tol * 100.0) {
                continue;
            }
            let cdist = bboxes[i].1.distance_to(&bboxes[j].1);
            out.push(FaceOverlap {
                face_a: fa,
                face_b: fb,
                same_surface: true,
                centroid_distance: cdist,
                surface_kind: surface_kind_label(&sa.geometry),
            });
        }
    }
    out
}

/// Collect the set of edges referenced by a face's outer loop.
fn face_edge_set(
    brep: &BRep,
    fh: FaceHandle,
) -> std::collections::BTreeSet<crate::topology::EdgeHandle> {
    let mut out = std::collections::BTreeSet::new();
    let Some(face) = brep.faces.get(fh) else {
        return out;
    };
    let Some(loop_data) = brep.loops.get(face.outer_loop) else {
        return out;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut safety = 256;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        out.insert(he.edge);
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    out
}

/// Detect half-edge pcurves that, when sampled, produce (u, v) values
/// outside the surface's valid parameter range. Slicers evaluating
/// such pcurves get garbage from the surface and treat the face as
/// non-manifold. Limited to surfaces with bounded UV (NURBS, sphere,
/// torus); planes and cylinders have unbounded ranges where escape
/// isn't well-defined.
fn detect_pcurve_escapes(brep: &BRep) -> Vec<PcurveEscape> {
    let mut out = Vec::new();
    let n_samples = 11;
    for (heh, he) in &brep.half_edges {
        let Some(pcurve) = he.pcurve.as_ref() else {
            continue;
        };
        let Some(loop_) = brep.loops.get(he.loop_) else {
            continue;
        };
        let face = loop_.face;
        let Some(surf) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let bounded = matches!(
            &surf.geometry,
            SurfaceGeometry::Nurbs(_)
                | SurfaceGeometry::Sphere { .. }
                | SurfaceGeometry::Torus { .. }
                | SurfaceGeometry::Gregory(_)
        );
        if !bounded {
            continue;
        }
        let ((u_lo, u_hi), (v_lo, v_hi)) = surf.geometry.parameter_range();
        if !u_lo.is_finite() || !u_hi.is_finite() || !v_lo.is_finite() || !v_hi.is_finite() {
            continue;
        }
        let (t_min, t_max) = pcurve.parameter_range();
        let mut max_overshoot = 0.0_f64;
        let mut worst_uv = Point3::new(0.0, 0.0, 0.0);
        let mut worst_t = t_min;
        for i in 0..n_samples {
            let s = i as f64 / (n_samples - 1) as f64;
            let t = t_min + s * (t_max - t_min);
            let uv = pcurve.point_at(t);
            let mut over = 0.0_f64;
            if uv.x < u_lo {
                over = over.max(u_lo - uv.x);
            }
            if uv.x > u_hi {
                over = over.max(uv.x - u_hi);
            }
            if uv.y < v_lo {
                over = over.max(v_lo - uv.y);
            }
            if uv.y > v_hi {
                over = over.max(uv.y - v_hi);
            }
            if over > max_overshoot {
                max_overshoot = over;
                worst_uv = uv;
                worst_t = t;
            }
        }
        if max_overshoot > 1e-9 {
            out.push(PcurveEscape {
                half_edge: heh,
                edge: he.edge,
                face,
                parameter: worst_t,
                uv: worst_uv,
                uv_range: ((u_lo, u_hi), (v_lo, v_hi)),
                max_overshoot,
            });
        }
    }
    out
}
