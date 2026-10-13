//! Hole-interior record-piece filter for the fragment walk (PK-428).
//!
//! Records deliberately run through a face's holes at full span — the
//! record clip is outer-loop-only, and hole semantics live in
//! containment. By walk time the graph has split those records at
//! junction nodes, and the pieces INSIDE one of the face's brought-in
//! holes are phantoms FOR THAT FACE: they bound none of its material,
//! yet both directions register to it and the chart walk chains them
//! into zero-area doubled rings (measured on the keyed-bore cap: a
//! 6-half-edge ring lying wholly inside the bore walked as an "inner
//! loop", the cap served with 3 inner rings, and cellular assembly
//! minted six degenerate junction stubs plus a once-used minor arc —
//! non-manifold, so the boolean fell to graph-merge, which walks the
//! same rings and served 74 open mesh edges).
//!
//! A Shared (record) half-edge whose midpoint lies strictly inside
//! one of the face's own holes is dropped from that face's walk list.
//! The face is the ORIGINAL operand face (the walk's breps predate
//! this boolean's splits), so its holes are exactly the brought-in
//! ones; the same graph half-edge stays registered to the other
//! operand's face, where it is a genuine ring member (a bore-interior
//! piece bounds the tool wall's in-bore fragment).
//!
//! Scope: faces carrying holes on planar charts (cap faces of bored
//! stock) and on singly-periodic analytic charts (a tube wall with a
//! drilled radial bore — measured serving 89 open mesh edges,
//! validator-clean, on the curved keyway sibling), with the query
//! branch-normalized per hole contour and seam-spanning contours
//! skipped. Doubly-periodic charts are excluded: their containment is
//! chart-cell territory, and a wrong fold would silently unwalk real
//! rings.
//!
//! The in-hole test samples each hole contour at 256 stations per
//! pcurve (the closed-bore precedent from the open-chain boundary
//! sampling): a 32-gon's sagitta (~5e-3 of a bore radius) would
//! swallow rim-hugging record pieces, while at 256 the polygon error
//! (~1e-5 of the radius) sits below the on-boundary band. Pieces
//! within the band of a hole rim are never filtered.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, Operand};

/// Stations per hole pcurve for the containment polygon.
const HOLE_CONTOUR_STATIONS: usize = 256;

/// Drop Shared record pieces lying strictly inside one of the face's
/// own holes. Returns `None` when the filter does not apply (no
/// surface, non-planar chart, or a face without holes) — the caller
/// then walks the unfiltered list.
pub(super) fn drop_hole_interior_record_pieces(
    graph: &IntersectionGraph,
    face_ref: &FaceRef,
    he_indices: &[usize],
    brep_a: &BRep,
    brep_b: &BRep,
    tolerance: &Tolerance,
) -> Option<Vec<usize>> {
    let brep = match face_ref.operand {
        Operand::A => brep_a,
        Operand::B => brep_b,
    };
    let face_data = brep.faces.get(face_ref.face)?;
    if face_data.inner_loops.is_empty() {
        return None;
    }
    let surface = face_data.surface.as_ref()?;
    // Planar charts and SINGLY-periodic analytic charts (a tube wall
    // with a drilled hole — measured torn on the curved keyway
    // sibling). The query point below is branch-normalized per hole
    // contour, and seam-spanning contours are skipped; doubly-periodic
    // charts stay excluded (their containment is chart-cell territory,
    // and a wrong fold here would silently unwalk real rings).
    let (u_period, v_period) = crate::clip::classification_chart_periods(&surface.geometry);
    let chart_ok = matches!(surface.geometry, SurfaceGeometry::Plane { .. })
        || (u_period.is_some() && v_period.is_none());
    if !chart_ok {
        return None;
    }
    // Hole contours as dense UV polygons; seam-spanning contours
    // (u-span over half the period — a hole the chart tears) are
    // skipped rather than mis-tested.
    let contours = crate::clip::face_uv_inner_boundaries(brep, face_ref.face);
    if contours.is_empty() {
        return None;
    }
    // Record midpoints first: a face no record piece reaches has
    // nothing to filter, and a hole whose padded boundary box holds no
    // midpoint cannot contain one, so its dense polygon is never built
    // (PK-802: a disc face with 63 slot holes sampled 64k stations per
    // cut to test the few pieces of the one new slot).
    let mids: Vec<Point3> = he_indices
        .iter()
        .filter_map(|&idx| record_midpoint(graph, idx))
        .collect();
    if mids.is_empty() {
        return None;
    }
    let periodic = u_period.is_some();
    let mut polygons: Vec<Vec<(f64, f64)>> = Vec::new();
    for contour in &contours {
        if !contour_box_holds_any(brep, contour, &mids, periodic, tolerance) {
            continue;
        }
        let mut polygon = Vec::new();
        for (_he, pcurve) in contour {
            let (t0, t1) = pcurve.parameter_range();
            if (t1 - t0).abs() < tolerance.numerical_floor() {
                continue;
            }
            for i in 0..HOLE_CONTOUR_STATIONS {
                let t = t0 + (t1 - t0) * i as f64 / HOLE_CONTOUR_STATIONS as f64;
                let q = pcurve.point_at(t);
                polygon.push((q.x, q.y));
            }
        }
        if polygon.len() < 3 {
            continue;
        }
        if let Some(p) = u_period {
            let (lo, hi) = polygon
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), q| {
                    (lo.min(q.0), hi.max(q.0))
                });
            if hi - lo > 0.5 * p {
                continue;
            }
        }
        polygons.push(polygon);
    }
    if polygons.is_empty() {
        return None;
    }
    let band = tolerance.linear() * 20.0;
    let strictly_in_hole = |u_raw: f64, v: f64| -> bool {
        for polygon in &polygons {
            // Branch-normalize the query to this contour: the record's
            // projection lands on the principal branch while the hole
            // pcurves may live a whole period away.
            let u = match u_period {
                Some(p) => {
                    let centre = polygon.iter().map(|q| q.0).sum::<f64>() / polygon.len() as f64;
                    u_raw + ((centre - u_raw) / p).round() * p
                }
                None => u_raw,
            };
            let mut inside = false;
            let mut min_d_sq = f64::INFINITY;
            let n = polygon.len();
            for i in 0..n {
                let (x0, y0) = polygon[i];
                let (x1, y1) = polygon[(i + 1) % n];
                if (y0 > v) != (y1 > v) {
                    let x_cross = x0 + (v - y0) / (y1 - y0) * (x1 - x0);
                    if x_cross > u {
                        inside = !inside;
                    }
                }
                let (dx, dy) = (x1 - x0, y1 - y0);
                let len_sq = dx * dx + dy * dy;
                let t = if len_sq > 1e-30 {
                    (((u - x0) * dx + (v - y0) * dy) / len_sq).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let (px, py) = (x0 + dx * t, y0 + dy * t);
                min_d_sq = min_d_sq.min((u - px) * (u - px) + (v - py) * (v - py));
            }
            if inside && min_d_sq > band * band {
                return true;
            }
        }
        false
    };
    let mut kept: Vec<usize> = Vec::with_capacity(he_indices.len());
    let mut dropped = 0usize;
    for &idx in he_indices {
        let he = &graph.half_edges[idx];
        if !matches!(he.source, EdgeSource::Shared { .. }) {
            kept.push(idx);
            continue;
        }
        let Some(mid) = record_midpoint(graph, idx) else {
            kept.push(idx);
            continue;
        };
        let Some((u, v)) =
            popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, &mid)
        else {
            kept.push(idx);
            continue;
        };
        if strictly_in_hole(u, v) {
            dropped += 1;
            if popkern_core::trace::trace_enabled("fragment-walk") {
                eprintln!(
                    "    [region-filter] {face_ref:?} drop he {idx} mid ({:.4},{:.4},{:.4}) inside hole",
                    mid.x, mid.y, mid.z
                );
            }
        } else {
            kept.push(idx);
        }
    }
    if dropped == 0 {
        return None;
    }
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [region-filter] {face_ref:?}: dropped {dropped} hole-interior record piece(s)"
        );
    }
    Some(kept)
}

/// Midpoint of a Shared (record) half-edge; `None` for any other.
///
/// A twin half-edge carries no curve of its own — its partner does,
/// and the curve midpoint is direction-independent. The chord fallback
/// is a last resort ONLY when neither direction has a curve: the chord
/// of a wide arc dips off the curve by the full sagitta, and on the nut
/// cap the tool tip-band's r = 3.15 arc (true midpoint azimuth −π/2,
/// far outside the bore) chord-sampled at r = 1.90 INSIDE the bore —
/// the twin dropped as hole-interior while its partner stayed, and the
/// 2-in/1-out imbalance declined both caps to the 3D-Newell fallback
/// (PK-464).
fn record_midpoint(graph: &IntersectionGraph, idx: usize) -> Option<Point3> {
    let he = &graph.half_edges[idx];
    if !matches!(he.source, EdgeSource::Shared { .. }) {
        return None;
    }
    let curve = he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref());
    Some(match curve {
        Some(c) => c.geometry.point_at(0.5 * (c.t_start + c.t_end)),
        None => {
            let a = graph.nodes[he.origin].position;
            let b = graph.nodes[graph.half_edges[he.twin].origin].position;
            Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
        }
    })
}

/// Whether any of `mids` lies in the 3D box of a hole contour's edge
/// curves, padded so the box holds the surface patch the contour
/// bounds: by the rim band on a plane (the patch is the contour's own
/// planar span), by the box diagonal on a periodic chart (the patch of
/// a contour spanning under half a period bulges off its boundary by
/// less than its chord). A contour with a curveless edge always holds.
fn contour_box_holds_any(
    brep: &BRep,
    contour: &[(
        popkern_core::topology::HalfEdgeHandle,
        popkern_core::geometry::pcurve::Pcurve,
    )],
    mids: &[Point3],
    periodic: bool,
    tolerance: &Tolerance,
) -> bool {
    let mut lo = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for (he, _) in contour {
        let Some(curve) = brep
            .half_edges
            .get(*he)
            .and_then(|h| brep.edges.get(h.edge))
            .and_then(|e| e.curve.as_ref())
        else {
            return true;
        };
        let b = curve.geometry.curve_bounding_box((
            curve.t_start.min(curve.t_end),
            curve.t_start.max(curve.t_end),
        ));
        lo = Point3::new(lo.x.min(b.min.x), lo.y.min(b.min.y), lo.z.min(b.min.z));
        hi = Point3::new(hi.x.max(b.max.x), hi.y.max(b.max.y), hi.z.max(b.max.z));
    }
    if !(lo.x.is_finite() && hi.x.is_finite()) {
        return true;
    }
    let diagonal = lo.distance_to(&hi);
    let pad = tolerance.linear() * 40.0 + if periodic { diagonal } else { 0.0 };
    mids.iter().any(|m| {
        m.x >= lo.x - pad
            && m.x <= hi.x + pad
            && m.y >= lo.y - pad
            && m.y <= hi.y + pad
            && m.z >= lo.z - pad
            && m.z <= hi.z + pad
    })
}
