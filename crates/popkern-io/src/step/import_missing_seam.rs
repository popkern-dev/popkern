//! Missing seams on closed NURBS bands (PK-641, PK-646).
//!
//! A face on a NURBS surface closed in one direction (its first and
//! last control rows coincide) is normally bounded through a SEAM
//! edge: the edge where the chart's two ends meet, used twice by the
//! face. Some writers bound a closed band by its two rims alone (an
//! extruded ellipse tube, a revolved collar: one loop per end, no
//! seam). The chart then has two rims that each run the full period
//! and nothing joining them — no charted polygon exists, and the face
//! meshed on the fallback grid, open along both rims (measured:
//! jack_in_the_box.stp, seven bodies).
//!
//! The repair adds the seam the topology lacks, by Euler operators:
//! each rim gets a vertex on the seam iso-curve (its own vertex when it
//! already lies there, otherwise the rim edge is split there), and a
//! make-edge-kill-ring joins the two, merging the rims into one loop.
//! The seam's curve is the surface's exact iso-curve at the chart end
//! (a straight ruling on an extrusion, a meridian section on a
//! revolve), trimmed between the two seam vertices.
//!
//! Pcurves: the rims were projected before the seam existed, so each
//! crossed the chart's wrap inside itself (a branch jump, or a closed
//! rim collapsed to a zero-span point); split at the seam, no rim
//! crosses the wrap any more, and the face's rim pcurves — and the
//! neighbouring faces' pcurves on split rim edges — are re-synthesized
//! from their curves. Each seam use runs along the chart end its
//! predecessor arrives at, so the loop is continuous in the chart.
//! The repair runs on a copy and commits only when the face
//! validates; every repair and decline is logged.

use popkern_core::geometry::curve::{
    CurveEval, CurveGeometry, CurveOnEdge, closest_point_on_curve,
};
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle,
};

use crate::step::diagnostics::{DiagnosticLog, Phase};

/// Add the seam to every closed-NURBS band face that lacks one;
/// returns the number of faces repaired.
pub(crate) fn add_missing_seams(brep: &mut BRep, log: &mut DiagnosticLog) -> usize {
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut repaired = 0;
    for face in faces {
        let Some(band) = band_surface(brep, face) else {
            continue;
        };
        // A closed rim's split direction is decided by a pcurve witness
        // inside `split_edge`; when it reads the rim the other way, the
        // split lands mirrored off the seam — retry mirrored.
        let mut candidate = brep.clone();
        let mut outcome = seam_face(&mut candidate, face, &band, false);
        if outcome.is_err() || candidate.validate().is_err() {
            let mut mirrored = brep.clone();
            if seam_face(&mut mirrored, face, &band, true).is_ok() && mirrored.validate().is_ok() {
                candidate = mirrored;
                outcome = Ok(());
            }
        }
        match outcome {
            Ok(()) => match candidate.validate() {
                Ok(()) => {
                    *brep = candidate;
                    repaired += 1;
                }
                Err(e) => log.warn(
                    Phase::MapTopology,
                    None,
                    format!(
                        "missing seam on {face:?} not added: the seamed face failed validation: {e}"
                    ),
                ),
            },
            Err(reason) => log.warn(
                Phase::MapTopology,
                None,
                format!("missing seam on {face:?} not added: {reason}"),
            ),
        }
    }
    if repaired > 0 {
        log.warn(
            Phase::MapTopology,
            None,
            format!("added the missing seam to {repaired} closed band face(s)"),
        );
    }
    repaired
}

/// A closed band's surface and the direction it closes in.
struct Band {
    nurbs: NurbsSurface,
    /// Closed in v (the seam is a u-running iso-curve) rather than u.
    closed_in_v: bool,
}

impl Band {
    /// The closed direction's parameter of a chart point.
    fn along(&self, uv: (f64, f64)) -> f64 {
        if self.closed_in_v { uv.1 } else { uv.0 }
    }

    /// The chart point at closed parameter `c`, across parameter `a`.
    fn chart(&self, c: f64, a: f64) -> (f64, f64) {
        if self.closed_in_v { (a, c) } else { (c, a) }
    }

    /// The closed direction's domain.
    fn period_range(&self) -> (f64, f64) {
        let ((u0, u1), (v0, v1)) = self.nurbs.parameter_range();
        if self.closed_in_v { (v0, v1) } else { (u0, u1) }
    }

    /// The exact iso-curve at the chart end, parameterized across.
    fn seam_curve(&self) -> NurbsCurve {
        let (c0, _) = self.period_range();
        if self.closed_in_v {
            self.nurbs.iso_curve_at_v(c0)
        } else {
            self.nurbs.iso_curve_at_u(c0)
        }
    }
}

/// The face's band when it needs a seam: a NURBS surface closed in
/// exactly one direction, two loops, and no edge used twice.
fn band_surface(brep: &BRep, face: FaceHandle) -> Option<Band> {
    let data = brep.faces.get(face)?;
    let SurfaceGeometry::Nurbs(nurbs) = &data.surface.as_ref()?.geometry else {
        return None;
    };
    if data.inner_loops.len() != 1 {
        return None;
    }
    let closed_in_v = match (closed_rows(nurbs, false), closed_rows(nurbs, true)) {
        (true, false) => false,
        (false, true) => true,
        _ => return None,
    };
    let mut edges = Vec::new();
    for lh in std::iter::once(data.outer_loop).chain(data.inner_loops.iter().copied()) {
        for he in loop_half_edges(brep, lh) {
            let edge = brep.half_edges.get(he)?.edge;
            if edges.contains(&edge) {
                return None; // already seamed
            }
            edges.push(edge);
        }
    }
    Some(Band {
        nurbs: nurbs.clone(),
        closed_in_v,
    })
}

/// The first and last control rows (in u, or in v) coincide.
fn closed_rows(n: &NurbsSurface, in_v: bool) -> bool {
    let scale = n
        .control_points
        .iter()
        .fold(0.0_f64, |m, p| {
            m.max(p.x.abs()).max(p.y.abs()).max(p.z.abs())
        })
        .max(1.0);
    let matches = |a: usize, b: usize| {
        n.control_points[a].distance_to(&n.control_points[b]) <= 1e-9 * scale
            && (n.weights[a] - n.weights[b]).abs() <= 1e-12 * n.weights[a].abs().max(1.0)
    };
    if in_v {
        n.num_v >= 2 && (0..n.num_u).all(|i| matches(i * n.num_v, i * n.num_v + n.num_v - 1))
    } else {
        let last = n.num_u - 1;
        n.num_u >= 2 && (0..n.num_v).all(|j| matches(j, last * n.num_v + j))
    }
}

fn loop_half_edges(brep: &BRep, lh: LoopHandle) -> Vec<HalfEdgeHandle> {
    let mut out = Vec::new();
    let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
        return out;
    };
    let mut he = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(he);
        match brep.half_edges.get(he) {
            Some(d) if d.next != start => he = d.next,
            _ => break,
        }
    }
    out
}

fn seam_face(brep: &mut BRep, face: FaceHandle, band: &Band, mirrored: bool) -> Result<(), String> {
    let surface = SurfaceGeometry::Nurbs(band.nurbs.clone());
    let (c0, c1) = band.period_range();
    let seam = band.seam_curve();
    // On the seam: within tolerance of the seam iso-curve (a
    // parametric test at the chart's end missed a rim vertex sitting
    // exactly on it — its projection is ambiguous between the ends).
    // The result is the across parameter.
    let gate = 10.0 * brep.tolerance.linear();
    let on_seam = |p: &Point3| -> Option<f64> {
        let t = closest_point_on_curve(&seam, p)?;
        (seam.point_at(t).distance_to(p) <= gate).then_some(t)
    };
    let data = brep.faces.get(face).ok_or("face vanished")?;
    let (outer, inner) = (data.outer_loop, data.inner_loops[0]);
    let (mut seam_he, mut split_edges) = (Vec::new(), Vec::new());
    for lh in [outer, inner] {
        let (he, split) = seam_vertex_on(brep, lh, &surface, band, &on_seam, mirrored)?;
        seam_he.push(he);
        split_edges.extend(split.into_iter().flatten());
    }
    let (va, vb) = (
        brep.half_edges[seam_he[0]].origin,
        brep.half_edges[seam_he[1]].origin,
    );
    let (a, b) = (brep.vertices[va].position, brep.vertices[vb].position);
    let (ta, tb) = (
        on_seam(&a).ok_or("outer seam vertex is off the seam iso-curve")?,
        on_seam(&b).ok_or("inner seam vertex is off the seam iso-curve")?,
    );
    let made = brep
        .mekr(seam_he[0], seam_he[1])
        .map_err(|e| format!("mekr: {e}"))?;
    brep.edges[made.edge].curve = Some(seam_geometry(brep, made.edge, &seam, (va, ta), (vb, tb)));
    // The seam's two uses hold placeholders while the rims' pcurves
    // are re-synthesized (measured on jack body 56: split rims kept
    // their pre-seam branch jumps, recorded as edge tolerances up to
    // 3.15; body 27's closed rims carried zero-span pcurves).
    let placeholder = Pcurve::linear(band.chart(c0, ta), band.chart(c0, tb), 0.0, 1.0);
    let seam_uses = [made.half_edge, made.half_edge_twin];
    for he in seam_uses {
        brep.half_edges[he].pcurve = Some(placeholder.clone());
    }
    for he in loop_half_edges(brep, brep.half_edges[made.half_edge].loop_) {
        if !seam_uses.contains(&he) {
            brep.half_edges[he].pcurve = None;
        }
    }
    brep.recompute_face_pcurves(face)
        .map_err(|e| format!("rim pcurves: {e}"))?;
    resynthesize_neighbours(brep, face, &split_edges)?;
    // Each seam use runs along the chart end its predecessor reaches.
    for (he, (from, to)) in [(made.half_edge, (ta, tb)), (made.half_edge_twin, (tb, ta))] {
        let prev = loop_half_edges(brep, brep.half_edges[he].loop_)
            .into_iter()
            .find(|&h| brep.half_edges[h].next == he)
            .ok_or("seam use has no predecessor")?;
        let pc = brep.half_edges[prev]
            .pcurve
            .as_ref()
            .ok_or("predecessor has no pcurve")?;
        let end = pc.point_at(pc.t_end);
        let reached = band.along((end.x, end.y));
        let side = if (reached - c0).abs() <= (reached - c1).abs() {
            c0
        } else {
            c1
        };
        brep.half_edges[he].pcurve = Some(Pcurve::linear(
            band.chart(side, from),
            band.chart(side, to),
            0.0,
            1.0,
        ));
    }
    Ok(())
}

/// The seam edge's curve, running from the edge's own half-edge's
/// origin. A straight iso-curve (an extrusion's ruling) is the exact
/// segment between the two vertices; a curved one (a revolve's
/// meridian section) is the iso-curve trimmed at the vertices'
/// parameters, its ends within the on-seam gate of the vertices, and
/// that end gap is recorded on the vertices and the edge (kernel
/// containment: vertex ≤ edge).
fn seam_geometry(
    brep: &mut BRep,
    edge: EdgeHandle,
    seam: &NurbsCurve,
    (va, ta): (VertexHandle, f64),
    (vb, tb): (VertexHandle, f64),
) -> CurveOnEdge {
    use popkern_core::precision::EntityHandle;
    let tol = brep.tolerance.linear();
    let first = brep.edges[edge].half_edge;
    let from_a = brep.half_edges.get(first).is_none_or(|h| h.origin == va);
    let ((vp, tp), (vq, tq)) = if from_a {
        ((va, ta), (vb, tb))
    } else {
        ((vb, tb), (va, ta))
    };
    let (p, q) = (brep.vertices[vp].position, brep.vertices[vq].position);
    if seam.degree == 1
        && seam.control_points.len() == 2
        && let Some(segment) = CurveOnEdge::line_segment(p, q, tol)
    {
        return segment;
    }
    let curve = if tp <= tq {
        CurveOnEdge {
            geometry: CurveGeometry::Nurbs(seam.clone()),
            t_start: tp,
            t_end: tq,
        }
    } else {
        // Reverse the iso-curve so its parameter runs p → q.
        let (k0, k1) = seam.parameter_range();
        let mut knots: Vec<f64> = seam.knots.iter().rev().map(|&k| k0 + k1 - k).collect();
        knots.sort_by(f64::total_cmp);
        let reversed = NurbsCurve {
            degree: seam.degree,
            knots,
            control_points: seam.control_points.iter().rev().copied().collect(),
            weights: seam.weights.iter().rev().copied().collect(),
        };
        CurveOnEdge {
            geometry: CurveGeometry::Nurbs(reversed),
            t_start: k0 + k1 - tp,
            t_end: k0 + k1 - tq,
        }
    };
    let gaps = [
        curve.geometry.point_at(curve.t_start).distance_to(&p),
        curve.geometry.point_at(curve.t_end).distance_to(&q),
    ];
    let mut widest = 0.0_f64;
    for (v, gap) in [(vp, gaps[0]), (vq, gaps[1])] {
        if gap > brep.vertex_tolerance(v) {
            brep.set_entity_tolerance(EntityHandle::Vertex(v), gap);
            widest = widest.max(gap);
        }
    }
    if widest > brep.edge_tolerance(edge) {
        brep.set_entity_tolerance(EntityHandle::Edge(edge), widest);
    }
    curve
}

/// A split rim's pieces carry, on the neighbouring face, halves of a
/// pcurve cut at the image of the split point; re-synthesize them
/// there, so both uses of every new edge are projections of its own
/// curve.
fn resynthesize_neighbours(
    brep: &mut BRep,
    face: FaceHandle,
    split_edges: &[EdgeHandle],
) -> Result<(), String> {
    let mut neighbours = Vec::new();
    for &edge in split_edges {
        let Some(first) = brep.edges.get(edge).map(|e| e.half_edge) else {
            continue;
        };
        for he in [first, brep.half_edges[first].twin] {
            let lh = brep.half_edges[he].loop_;
            let Some(owner) = brep.loops.get(lh).map(|l| l.face) else {
                continue;
            };
            if owner != face {
                brep.half_edges[he].pcurve = None;
                neighbours.push(owner);
            }
        }
    }
    for owner in neighbours {
        brep.recompute_face_pcurves(owner)
            .map_err(|e| format!("neighbour pcurves: {e}"))?;
    }
    Ok(())
}

/// A half-edge of loop `lh` whose origin lies on the seam iso-curve,
/// splitting the loop's rim edge at its seam crossing when no vertex
/// does (the split's two edges are returned).
fn seam_vertex_on(
    brep: &mut BRep,
    lh: LoopHandle,
    surface: &SurfaceGeometry,
    band: &Band,
    on_seam: &dyn Fn(&Point3) -> Option<f64>,
    mirrored: bool,
) -> Result<(HalfEdgeHandle, Option<[EdgeHandle; 2]>), String> {
    let half_edges = loop_half_edges(brep, lh);
    if let Some(&he) = half_edges
        .iter()
        .find(|&&he| on_seam(&brep.vertices[brep.half_edges[he].origin].position).is_some())
    {
        return Ok((he, None));
    }
    let (c0, c1) = band.period_range();
    let span = c1 - c0;
    for he in half_edges {
        let edge = brep.half_edges[he].edge;
        let Some(curve) = brep.edges[edge].curve.clone() else {
            continue;
        };
        let along_at = |t: f64| {
            surface
                .closest_point(&curve.geometry.point_at(t))
                .map(|uv| band.along(uv))
        };
        let n = 64;
        let ts: Vec<f64> = (0..=n)
            .map(|k| curve.t_start + (curve.t_end - curve.t_start) * k as f64 / n as f64)
            .collect();
        for w in ts.windows(2) {
            let (Some(ca), Some(cb)) = (along_at(w[0]), along_at(w[1])) else {
                continue;
            };
            if (ca - cb).abs() < 0.5 * span {
                continue; // no wrap across this sample interval
            }
            // Bisect toward the chart wrap: the side flips at the seam.
            let low_side = |c: f64| c - c0 < 0.5 * span;
            let (mut lo, mut hi) = (w[0], w[1]);
            let side_lo = low_side(ca);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                match along_at(mid) {
                    Some(c) if low_side(c) == side_lo => lo = mid,
                    Some(_) => hi = mid,
                    None => break,
                }
            }
            // Split at the bisected curve parameter itself (re-inverting
            // the point on a closed ellipse rim can fail).
            let fraction = (0.5 * (lo + hi) - curve.t_start) / (curve.t_end - curve.t_start);
            let fraction = if mirrored { 1.0 - fraction } else { fraction };
            let split = brep
                .split_edge(edge, fraction)
                .map_err(|e| format!("rim split: {e}"))?;
            return loop_half_edges(brep, lh)
                .into_iter()
                .find(|&h| brep.half_edges[h].origin == split.vertex)
                .map(|h| (h, Some([split.edge_a, split.edge_b])))
                .ok_or_else(|| "split vertex not on its loop".to_string());
        }
    }
    Err("no rim crossing of the seam found".into())
}
