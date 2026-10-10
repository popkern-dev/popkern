//! Tangent-chain junction sections in per-edge fillets (PK-770).
//!
//! Where two fillet strips meet at a tangent-chain junction, the
//! per-edge pipeline anchors both strips on the junction's exact
//! section contacts (`g1_junction::junction_section_contact`). The
//! corner vertex between them is moved onto the section arc by the
//! terminal settlement, but the two stub edges joining it to the
//! contacts stay straight chords. Measured on a rounded rim's
//! line → arc junction at r = 3: 0.228 off both strips, exactly the
//! sagitta of a 45° chord. Nothing else authored the section: the
//! closed and open chain surgery builds its own; this is the per-edge
//! analogue.
//!
//! The same signature appears on a strip's RAIL split at a valence-2
//! vertex (a torus strip's contact with a corner cylinder, split at the
//! cylinder's seam): its fitted pieces sat 1.7e-3 off, and the circle
//! through their three points is the exact contact circle.
//!
//! **Repair, by construction.** A vertex M on exactly two edges, both
//! between the same two CREATED faces (strips, or supports the blend
//! split), is a junction corner or a split rail. Its two stubs' far ends and M lie on the section circle, and
//! three points determine it. The circle is accepted only when its arcs
//! lie on BOTH strips: the arc midpoints are within the acceptance band
//! of each strip surface. Both stubs are then re-authored as the exact
//! short arcs and the two strips' pcurves recomputed. Zero topology
//! change; any doubt leaves the stubs as they were (the curve-surface
//! gate still reports them).
//!
//! **Tolerance:** 100 × linear for the on-strip acceptance (the same
//! band as the chain surgery's), local to this pass.

use std::collections::HashSet;

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// The circle through three points: centre, unit normal, radius.
fn circle_through(a: Point3, b: Point3, c: Point3) -> Option<(Point3, Vec3, f64)> {
    let (ab, ac) = (b - a, c - a);
    let n = ab.cross(&ac);
    let n2 = n.dot(&n);
    // Nearly collinear points define no section (measured: a straight
    // rail's three points gave an 8e15-radius "arc").
    if n2 <= (1e-6 * ab.length() * ac.length()).powi(2) {
        return None;
    }
    let to_center = (n.cross(&ab) * ac.dot(&ac) + ac.cross(&n) * ab.dot(&ab)) * (0.5 / n2);
    let center = a + to_center;
    Some((center, n.normalized()?, to_center.length()))
}

/// The two faces of an edge.
fn faces_of(brep: &BRep, edge: EdgeHandle) -> Option<(FaceHandle, FaceHandle)> {
    let he = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let tw = brep.half_edges.get(he.twin)?;
    Some((
        brep.loops.get(he.loop_)?.face,
        brep.loops.get(tw.loop_)?.face,
    ))
}

/// The edge's other endpoint.
fn far_end(brep: &BRep, edge: EdgeHandle, from: VertexHandle) -> Option<VertexHandle> {
    let he = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let to = brep.half_edges.get(he.twin)?.origin;
    Some(if he.origin == from { to } else { he.origin })
}

/// Settle every junction corner between two strips the blend created
/// (faces of `result` absent from `source`). Returns the number of
/// junction corners whose stubs were re-authored as section arcs.
pub(super) fn settle_junction_sections(
    source: &BRep,
    result: &mut BRep,
    tolerance: &Tolerance,
) -> usize {
    let created: HashSet<FaceHandle> = result
        .faces
        .keys()
        .filter(|f| !source.faces.contains_key(*f))
        .collect();
    let band = tolerance.linear().max(1e-9) * 100.0;
    // Vertex → its edges.
    let mut at: std::collections::BTreeMap<VertexHandle, Vec<EdgeHandle>> = Default::default();
    for (e, data) in &result.edges {
        let Some(he) = result.half_edges.get(data.half_edge) else {
            continue;
        };
        let Some(tw) = result.half_edges.get(he.twin) else {
            continue;
        };
        at.entry(he.origin).or_default().push(e);
        if tw.origin != he.origin {
            at.entry(tw.origin).or_default().push(e);
        }
    }
    let mut settled = 0;
    for (m, edges) in at {
        let [e1, e2] = edges.as_slice() else {
            continue;
        };
        let (Some(f1), Some(f2)) = (faces_of(result, *e1), faces_of(result, *e2)) else {
            continue;
        };
        let strips = |(a, b): (FaceHandle, FaceHandle)| {
            a != b && created.contains(&a) && created.contains(&b)
        };
        let same = f1 == f2 || f1 == (f2.1, f2.0);
        if !same || !strips(f1) {
            continue;
        }
        if let Some(n) = settle_corner(result, m, [*e1, *e2], f1, band) {
            settled += n;
        }
    }
    if settled > 0 && popkern_core::trace::trace_enabled("terminal-arc") {
        eprintln!("[junction-section] {settled} junction corner(s) settled on section arcs");
    }
    settled
}

/// One junction corner M with stubs `edges` between `strips`.
fn settle_corner(
    brep: &mut BRep,
    m: VertexHandle,
    edges: [EdgeHandle; 2],
    strips: (FaceHandle, FaceHandle),
    band: f64,
) -> Option<usize> {
    let pm = brep.vertices.get(m)?.position;
    let ends = [far_end(brep, edges[0], m)?, far_end(brep, edges[1], m)?];
    let (pa, pb) = (
        brep.vertices.get(ends[0])?.position,
        brep.vertices.get(ends[1])?.position,
    );
    let (center, normal, radius) = circle_through(pa, pm, pb)?;
    // Already exact: both stubs lie on both strips.
    let surface = |f: FaceHandle| Some(brep.faces.get(f)?.surface.as_ref()?.geometry.clone());
    let (s1, s2) = (surface(strips.0)?, surface(strips.1)?);
    let off = |s: &popkern_core::geometry::surface::SurfaceGeometry, p: Point3| {
        s.closest_point(&p)
            .map(|(u, v)| s.point_at(u, v).distance_to(&p))
            .unwrap_or(f64::INFINITY)
    };
    let mut arcs = Vec::with_capacity(2);
    for (edge, end) in edges.iter().zip([pa, pb]) {
        let arc = super::tangent_chain_geometry::section_arc(end, pm, center, normal, radius)?;
        let mid = arc.geometry.point_at(0.5 * (arc.t_start + arc.t_end));
        if off(&s1, mid) > band || off(&s2, mid) > band {
            return None;
        }
        let current = brep.edges.get(*edge)?.curve.as_ref()?;
        let cur_mid = current
            .geometry
            .point_at(0.5 * (current.t_start + current.t_end));
        if off(&s1, cur_mid) <= band && off(&s2, cur_mid) <= band {
            return None; // already on both strips: nothing to settle
        }
        arcs.push((*edge, arc));
    }
    for (edge, arc) in arcs {
        let data = brep.edges.get_mut(edge)?;
        data.curve = Some(arc);
        let he = data.half_edge;
        let tw = brep.half_edges.get(he)?.twin;
        for h in [he, tw] {
            if let Some(d) = brep.half_edges.get_mut(h) {
                d.pcurve = None;
            }
        }
    }
    // The edge's canonical half-edge must run with its curve: the arc
    // runs from the far end to M.
    for (edge, end) in edges.iter().zip(ends) {
        let he = brep.edges.get(*edge)?.half_edge;
        if brep.half_edges.get(he)?.origin != end {
            let tw = brep.half_edges.get(he)?.twin;
            brep.edges.get_mut(*edge)?.half_edge = tw;
        }
    }
    for f in [strips.0, strips.1] {
        brep.recompute_face_pcurves(f).ok()?;
    }
    Some(1)
}
