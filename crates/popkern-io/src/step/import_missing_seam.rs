//! Missing seams on closed NURBS bands (PK-641).
//!
//! A face on a NURBS surface closed in one direction (its first and
//! last control rows coincide) is normally bounded through a SEAM
//! edge: the edge where the chart's two ends meet, used twice by the
//! face. Some writers bound a closed band by its two rims alone (an
//! extruded ellipse tube: one loop per end, no seam). The chart then
//! has two rims that each run the full period and nothing joining
//! them — no charted polygon exists, and the face meshed on the
//! fallback grid, open along both rims (measured: jack_in_the_box.stp,
//! six bodies).
//!
//! The repair adds the seam the topology lacks, by Euler operators:
//! each rim gets a vertex on the seam iso-line (its own vertex when it
//! already lies there, otherwise the rim edge is split there), and a
//! make-edge-kill-ring joins the two, merging the rims into one loop.
//! The seam curve is the iso-line, a straight segment on a surface
//! ruled across the band (an extrusion); a curved iso-line declines.
//! Each seam half-edge's pcurve runs along the chart end its
//! predecessor arrives at, so the loop is continuous in the chart (the
//! rims keep their own pcurves, split with them). The
//! repair runs on a copy and commits only when the face validates;
//! every repair and decline is logged.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use crate::step::diagnostics::{DiagnosticLog, Phase};

/// Add the seam to every closed-NURBS band face that lacks one;
/// returns the number of faces repaired.
pub(crate) fn add_missing_seams(brep: &mut BRep, log: &mut DiagnosticLog) -> usize {
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut repaired = 0;
    for face in faces {
        let Some(nurbs) = band_surface(brep, face) else {
            continue;
        };
        // A closed rim's split direction is decided by a pcurve witness
        // inside `split_edge`; when it reads the rim the other way, the
        // split lands mirrored off the seam — retry mirrored.
        let mut candidate = brep.clone();
        let mut outcome = seam_face(&mut candidate, face, &nurbs, false);
        if outcome.is_err() || candidate.validate().is_err() {
            let mut mirrored = brep.clone();
            if seam_face(&mut mirrored, face, &nurbs, true).is_ok() && mirrored.validate().is_ok() {
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

/// The face's NURBS surface when it is a band in need of a seam:
/// closed in u, ruled (straight iso-lines) across, two loops, and no
/// edge used twice.
fn band_surface(brep: &BRep, face: FaceHandle) -> Option<NurbsSurface> {
    let data = brep.faces.get(face)?;
    let SurfaceGeometry::Nurbs(nurbs) = &data.surface.as_ref()?.geometry else {
        return None;
    };
    if data.inner_loops.len() != 1 || !closed_in_u(nurbs) || nurbs.degree_v != 1 || nurbs.num_v != 2
    {
        return None;
    }
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
    Some(nurbs.clone())
}

fn closed_in_u(n: &NurbsSurface) -> bool {
    let scale = n
        .control_points
        .iter()
        .fold(0.0_f64, |m, p| {
            m.max(p.x.abs()).max(p.y.abs()).max(p.z.abs())
        })
        .max(1.0);
    let last = n.num_u - 1;
    n.num_u >= 2
        && (0..n.num_v).all(|j| {
            let (a, b) = (j, last * n.num_v + j);
            n.control_points[a].distance_to(&n.control_points[b]) <= 1e-9 * scale
                && (n.weights[a] - n.weights[b]).abs() <= 1e-12 * n.weights[a].abs().max(1.0)
        })
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

fn seam_face(
    brep: &mut BRep,
    face: FaceHandle,
    nurbs: &NurbsSurface,
    mirrored: bool,
) -> Result<(), String> {
    let surface = SurfaceGeometry::Nurbs(nurbs.clone());
    let ((u0, u1), _) = nurbs.parameter_range();
    // On the seam: within tolerance of the seam iso-line, which on a
    // band ruled across is the segment between its column's two
    // control points (a parametric test at the chart's end missed a
    // rim vertex sitting exactly on it — its projection is ambiguous
    // between u0 and u1). The v is the projection's.
    let (p_lo, p_hi) = (nurbs.control_points[0], nurbs.control_points[1]);
    let gate = 10.0 * brep.tolerance.linear();
    let on_seam = |p: &Point3| -> Option<f64> {
        let d = p_hi - p_lo;
        let t = ((*p - p_lo).dot(&d) / d.dot(&d).max(1e-300)).clamp(0.0, 1.0);
        if p.distance_to(&(p_lo + d * t)) > gate {
            return None;
        }
        surface.closest_point(p).map(|(_, v)| v)
    };
    let data = brep.faces.get(face).ok_or("face vanished")?;
    let (outer, inner) = (data.outer_loop, data.inner_loops[0]);
    let mut seam_he = Vec::new();
    for lh in [outer, inner] {
        seam_he.push(seam_vertex_on(
            brep,
            lh,
            &surface,
            (u0, u1),
            &on_seam,
            mirrored,
        )?);
    }
    let a = brep.vertices[brep.half_edges[seam_he[0]].origin].position;
    let b = brep.vertices[brep.half_edges[seam_he[1]].origin].position;
    // The iso-line must be the straight segment between the two.
    let (va, vb) = (
        on_seam(&a).ok_or("outer seam vertex left the seam")?,
        on_seam(&b).ok_or("inner seam vertex left the seam")?,
    );
    // The seam iso-line u = u0 of a band ruled across (degree 1, two
    // control columns in v) is the segment between its two control
    // points; both seam vertices must lie on it.
    let (p_lo, p_hi) = (nurbs.control_points[0], nurbs.control_points[1]);
    let off_line = |p: &Point3| {
        let d = p_hi - p_lo;
        let t = (*p - p_lo).dot(&d) / d.dot(&d).max(1e-300);
        p.distance_to(&(p_lo + d * t))
    };
    let gate = 10.0 * brep.tolerance.linear();
    if off_line(&a) > gate || off_line(&b) > gate {
        return Err("the seam vertices are off the seam iso-line".into());
    }
    let made = brep
        .mekr(seam_he[0], seam_he[1])
        .map_err(|e| format!("mekr: {e}"))?;
    let tol = brep.tolerance.linear();
    brep.edges[made.edge].curve = CurveOnEdge::line_segment(a, b, tol);
    // The rims keep their pcurves (split with them); only the seam's
    // two uses are new. Recomputing the face's pcurves re-trimmed a
    // split closed rim onto the wrong 0/2π branch (a 3.4 shell gap).
    // Each seam use runs along the chart end its predecessor reaches.
    for (he, (from, to)) in [(made.half_edge, (va, vb)), (made.half_edge_twin, (vb, va))] {
        let prev = loop_half_edges(brep, brep.half_edges[he].loop_)
            .into_iter()
            .find(|&h| brep.half_edges[h].next == he)
            .ok_or("seam use has no predecessor")?;
        let pc = brep.half_edges[prev]
            .pcurve
            .as_ref()
            .ok_or("predecessor has no pcurve")?;
        let end = pc.point_at(pc.t_end);
        let u_side = if (end.x - u0).abs() <= (end.x - u1).abs() {
            u0
        } else {
            u1
        };
        brep.half_edges[he].pcurve = Some(Pcurve::linear((u_side, from), (u_side, to), 0.0, 1.0));
    }
    Ok(())
}

/// A half-edge of loop `lh` whose origin lies on the seam iso-line,
/// splitting the loop's rim edge at its seam crossing when no vertex
/// does.
fn seam_vertex_on(
    brep: &mut BRep,
    lh: LoopHandle,
    surface: &SurfaceGeometry,
    (u0, u1): (f64, f64),
    on_seam: &dyn Fn(&Point3) -> Option<f64>,
    mirrored: bool,
) -> Result<HalfEdgeHandle, String> {
    let half_edges = loop_half_edges(brep, lh);
    if let Some(&he) = half_edges
        .iter()
        .find(|&&he| on_seam(&brep.vertices[brep.half_edges[he].origin].position).is_some())
    {
        return Ok(he);
    }
    let span = u1 - u0;
    for he in half_edges {
        let edge = brep.half_edges[he].edge;
        let Some(curve) = brep.edges[edge].curve.clone() else {
            continue;
        };
        let u_at = |t: f64| {
            surface
                .closest_point(&curve.geometry.point_at(t))
                .map(|(u, _)| u)
        };
        let n = 64;
        let ts: Vec<f64> = (0..=n)
            .map(|k| curve.t_start + (curve.t_end - curve.t_start) * k as f64 / n as f64)
            .collect();
        for w in ts.windows(2) {
            let (Some(ua), Some(ub)) = (u_at(w[0]), u_at(w[1])) else {
                continue;
            };
            if (ua - ub).abs() < 0.5 * span {
                continue; // no wrap across this sample interval
            }
            // Bisect toward the chart wrap: the side flips at the seam.
            let low_side = |u: f64| u - u0 < 0.5 * span;
            let (mut lo, mut hi) = (w[0], w[1]);
            let side_lo = low_side(ua);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                match u_at(mid) {
                    Some(u) if low_side(u) == side_lo => lo = mid,
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
                .ok_or_else(|| "split vertex not on its loop".to_string());
        }
    }
    Err("no rim crossing of the seam found".into())
}
