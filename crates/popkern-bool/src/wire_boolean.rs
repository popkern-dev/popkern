//! Booleans between a wire and a solid (plan 4.12 item 1, PK-1068).
//!
//! A wire body is edges on no face (`BodyKind::Wire`). Against a solid
//! it is cut where its curves pierce or leave the solid's faces, and
//! each piece is kept by where it lies:
//!
//! - wire ∩ solid keeps the pieces inside the closed solid, those lying
//!   on its boundary included (the closure of the intersection);
//! - wire − solid keeps the pieces strictly outside;
//! - solid − wire is the solid itself: a wire bounds no volume, so
//!   regularized subtraction removes nothing;
//! - solid ∪ wire is the solid with the wire's outside pieces beside
//!   it, one mixed body (`BodyKind::Mixed`).
//!
//! The cut parameters are the curve's crossings with every face of the
//! solid's face surfaces (`intersect_surface_curve_dense`, the solver
//! boolean edge-hit discovery uses, so a multi-wrap face misses no
//! crossing), the ends of every stretch the curve lies on one, and its
//! meetings with the solid's edges. A cut outside a face's trim only
//! splits a piece whose halves lie alike; they join again. Each piece between consecutive cuts takes the
//! classification of its midpoint, which the cuts make uniform along
//! the piece: on the boundary when it lies on a face within the
//! tolerance (inside the trim), else `point_in_solid_final`'s verdict. The kept pieces keep the wire's
//! curves (trimmed ranges of the same geometry) and its connectivity: a
//! piece reaching an original vertex reuses one vertex for it, so two
//! wire edges kept up to their shared vertex stay joined. Cuts closer
//! than the linear tolerance in space are one cut.

use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::intersection::surface_curve::{
    SurfaceCurveIntersection, intersect_surface_curve_dense,
};
use popkern_core::geometry::intersection::{CurveIntersection, intersect_curves};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::bounded_parameter_range;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind, VertexHandle};

use crate::error::{BoolError, BoolOperation, PointClassification};

/// Where a piece of wire lies against the solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Lies {
    Inside,
    On,
    Outside,
}

/// A boolean between a wire body and a solid, in either order (module
/// doc for what each operation keeps).
///
/// # Errors
///
/// [`BoolError::UnsupportedBodyKind`] unless one operand is a wire and
/// the other a solid; [`BoolError::MergeFailed`] when a crossing or a
/// classification declines.
pub fn boolean_with_wire(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tolerance: Tolerance,
) -> Result<BRep, BoolError> {
    const ACCEPTED: &[BodyKind] = &[BodyKind::Wire, BodyKind::Solid];
    let kinds = (a.body_kind(), b.body_kind());
    let unsupported = |operand: usize, kind: BodyKind| BoolError::UnsupportedBodyKind {
        operation: "wire boolean",
        operand,
        kind,
        accepted: ACCEPTED,
    };
    let (wire, solid, wire_first) = match kinds {
        (BodyKind::Wire, BodyKind::Solid) => (a, b, true),
        (BodyKind::Solid, BodyKind::Wire) => (b, a, false),
        (BodyKind::Wire, other) => return Err(unsupported(1, other)),
        (other, _) => return Err(unsupported(0, other)),
    };
    match (op, wire_first) {
        (BoolOperation::Intersection, _) => {
            wire_part(wire, solid, tolerance, |l| l != Lies::Outside)
        }
        (BoolOperation::Subtraction, true) => {
            wire_part(wire, solid, tolerance, |l| l == Lies::Outside)
        }
        (BoolOperation::Subtraction, false) => Ok(solid.clone()),
        (BoolOperation::Union, _) => {
            let outside = wire_part(wire, solid, tolerance, |l| l == Lies::Outside)?;
            crate::operation::helpers::combine_breps(solid.clone(), outside)
        }
    }
}

fn declined(reason: String) -> BoolError {
    BoolError::MergeFailed {
        stage: "Classification".to_string(),
        reason: format!("wire boolean: {reason}"),
    }
}

/// The parameters in `(lo, hi)` where `curve` crosses or starts or ends
/// lying on a face of `solid`.
fn cuts(
    curve: &CurveOnEdge,
    (lo, hi): (f64, f64),
    solid: &BRep,
    tolerance: &Tolerance,
) -> Result<Vec<f64>, BoolError> {
    let mut out = Vec::new();
    for (fh, face) in &solid.faces {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        let positions = solid.face_outer_loop_positions(fh);
        let (u_range, v_range) = bounded_parameter_range(&positions, &surface.geometry);
        let hits = intersect_surface_curve_dense(
            &surface.geometry,
            u_range,
            v_range,
            &curve.geometry,
            (lo, hi),
            tolerance,
        )
        .map_err(|e| declined(format!("curve against face {fh:?}: {e}")))?;
        for hit in hits.intersections {
            match hit {
                // Every crossing of the face's surface cuts, trimmed or
                // not: a cut where the surface runs past the face splits
                // a piece whose two halves lie alike, and they join
                // again (a trim test rejected a ball's real crossings at
                // its seam: measured, a circle through it kept whole).
                SurfaceCurveIntersection::Point { t, .. } => out.push(t),
                SurfaceCurveIntersection::Overlap { t_range } => {
                    out.extend([t_range.0, t_range.1]);
                }
            }
        }
    }
    // Where the wire meets an edge of the solid: a stretch on a face
    // starts and ends on the face's boundary, where the face crossing
    // can fall a rounding outside the face's trim (measured: a wire
    // along a rotated cube's top-face edge kept no cut there).
    for (_, edge) in &solid.edges {
        let Some(e) = edge.curve.as_ref() else {
            continue;
        };
        let range = (e.t_start.min(e.t_end), e.t_start.max(e.t_end));
        if range.0 >= range.1 {
            continue;
        }
        let Ok(hits) = intersect_curves(&curve.geometry, (lo, hi), &e.geometry, range, tolerance)
        else {
            continue;
        };
        for hit in hits.intersections {
            match hit {
                CurveIntersection::Point { t1, .. } => out.push(t1),
                CurveIntersection::Overlap { t1_range, .. } => out.extend([t1_range.0, t1_range.1]),
            }
        }
    }
    Ok(out)
}

/// The pieces of `wire` whose lie `keep` accepts, as a wire body.
fn wire_part(
    wire: &BRep,
    solid: &BRep,
    tolerance: Tolerance,
    keep: impl Fn(Lies) -> bool,
) -> Result<BRep, BoolError> {
    let boxes = crate::intersect::compute_face_aabbs(solid);
    let linear = tolerance.linear();
    // On a face: within the tolerance of its surface, inside its trim.
    // Asked first: the ray classifier calls a point on a face outside
    // (measured: a wire along a box's top face).
    let on_boundary = |p: &popkern_core::geometry::Point3| {
        solid.faces.iter().any(|(fh, face)| {
            let Some(surface) = face.surface.as_ref() else {
                return false;
            };
            // The chart the face's trims are in (a plane's own
            // closest-point parameters differ once the body is rotated:
            // measured, a wire on a rotated cube's face read outside).
            let Some((u, v)) =
                popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, p)
            else {
                return false;
            };
            surface.geometry.point_at(u, v).distance_to(p) <= linear
                && crate::classify_uv_loops::point_in_face_region_exact(solid, fh, u, v)
                    .unwrap_or(false)
        })
    };
    let lies = |t: f64, curve: &CurveOnEdge| -> Result<Lies, BoolError> {
        let p = curve.geometry.point_at(t);
        if on_boundary(&p) {
            return Ok(Lies::On);
        }
        match crate::classify::point_in_solid_final(&p, solid, &boxes, &tolerance) {
            Ok(PointClassification::Inside) => Ok(Lies::Inside),
            Ok(PointClassification::Outside) => Ok(Lies::Outside),
            Ok(PointClassification::OnBoundary) => Ok(Lies::On),
            Err(e) => Err(declined(format!("classifying {p:?}: {e}"))),
        }
    };
    let mut builder = BRepBuilder::new(tolerance);
    // Output vertices standing for the wire's own vertices.
    let mut at_vertex: std::collections::HashMap<VertexHandle, VertexHandle> =
        std::collections::HashMap::new();
    for (_, edge) in &wire.edges {
        let Some(curve) = edge.curve.as_ref() else {
            return Err(declined("a wire edge has no curve".to_string()));
        };
        let he = &wire.half_edges[edge.half_edge];
        let start_vertex = he.origin;
        let end_vertex = he
            .mate()
            .and_then(|m| wire.half_edges.get(m))
            .map_or(start_vertex, |t| t.origin);
        let (lo, hi) = (
            curve.t_start.min(curve.t_end),
            curve.t_start.max(curve.t_end),
        );
        let mut interior: Vec<f64> = cuts(curve, (lo, hi), solid, &tolerance)?
            .into_iter()
            .filter(|t| *t > lo && *t < hi)
            .collect();
        interior.sort_by(f64::total_cmp);
        // One cut per place in space; the curve's own ends stay apart
        // even where they meet (a closed edge).
        let at = |t: f64| curve.geometry.point_at(t);
        let mut merged: Vec<f64> = vec![lo];
        for t in interior {
            if merged
                .last()
                .is_some_and(|&q| at(q).distance_to(&at(t)) <= linear)
            {
                continue;
            }
            merged.push(t);
        }
        let n = merged.len();
        if n > 1 && at(merged[n - 1]).distance_to(&at(hi)) <= linear {
            merged[n - 1] = hi;
        } else {
            merged.push(hi);
        }
        // Consecutive kept pieces join into one stretch.
        let mut stretches: Vec<(f64, f64)> = Vec::new();
        for w in merged.windows(2) {
            let (t0, t1) = (w[0], w[1]);
            if !keep(lies(0.5 * (t0 + t1), curve)?) {
                continue;
            }
            match stretches.last_mut() {
                Some(last) if last.1 == t0 => last.1 = t1,
                _ => stretches.push((t0, t1)),
            }
        }
        let forward = curve.t_start <= curve.t_end;
        // The wire vertex at parameter `t`, if it is an end of the edge.
        let end_at = |t: f64| {
            if t == curve.t_start {
                Some(start_vertex)
            } else if t == curve.t_end {
                Some(end_vertex)
            } else {
                None
            }
        };
        for (t0, t1) in stretches {
            let (from, to) = if forward { (t0, t1) } else { (t1, t0) };
            let mut vertex = |t: f64, builder: &mut BRepBuilder| match end_at(t) {
                Some(v) => *at_vertex
                    .entry(v)
                    .or_insert_with(|| builder.add_vertex(wire.vertices[v].position)),
                None => builder.add_vertex(curve.geometry.point_at(t)),
            };
            let v_from = vertex(from, &mut builder);
            // A whole closed edge keeps its one vertex.
            let closed = end_at(from).is_some() && end_at(from) == end_at(to);
            let v_to = if closed {
                v_from
            } else {
                vertex(to, &mut builder)
            };
            let (eh, _, _) = builder.add_edge(v_from, v_to);
            builder.set_edge_curve(
                eh,
                CurveOnEdge {
                    geometry: curve.geometry.clone(),
                    t_start: from,
                    t_end: to,
                },
            );
        }
    }
    Ok(builder.build_unchecked())
}

#[cfg(test)]
#[path = "wire_boolean_tests.rs"]
mod tests;
