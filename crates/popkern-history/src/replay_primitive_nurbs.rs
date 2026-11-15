//! The B-spline sheet primitive (`PrimitiveShape::SheetNurbs`): the face
//! of a surface given by its control net, over its whole chart.
//!
//! The net, weights and knots are checked as a B-spline surface needs
//! them (a rectangular net of at least `degree + 1` points each way,
//! positive weights, clamped non-decreasing knots of the matching
//! count), the surface is built as authored, and its sheet comes from
//! the sweep topology builder for an open profile: one face, its four
//! boundary iso-curves free edges. A net that closes on itself (a row
//! collapsed to a pole, two rows coinciding) has no free boundary there
//! and declines.

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{failed, invalid};

#[path = "replay_primitive_nurbs_ball.rs"]
mod replay_primitive_nurbs_ball;
use crate::error::ReplayError;

/// The sheet of the B-spline surface with control net `points`.
pub(super) fn sheet(
    shape: &str,
    points: &[Vec<[f64; 3]>],
    weights: &[Vec<f64>],
    degree: [u32; 2],
    knots: &[Vec<f64>; 2],
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let surface = surface(shape, points, weights, degree, knots)?;
    // An edge of the chart collapsed to a point, or two opposite edges
    // coinciding, closes the sheet there: no free boundary.
    let (u0, u1) = (surface.knots_u[0], *surface.knots_u.last().unwrap_or(&1.0));
    let (v0, v1) = (surface.knots_v[0], *surface.knots_v.last().unwrap_or(&1.0));
    use popkern_core::geometry::surface::SurfaceEval;
    let iso_spread = |fixed_u: Option<f64>, fixed_v: Option<f64>| -> f64 {
        let pts: Vec<Point3> = (0..=8)
            .map(|k| {
                let t = k as f64 / 8.0;
                match (fixed_u, fixed_v) {
                    (Some(u), _) => surface.point_at(u, v0 + t * (v1 - v0)),
                    (_, Some(v)) => surface.point_at(u0 + t * (u1 - u0), v),
                    _ => surface.point_at(u0, v0),
                }
            })
            .collect();
        pts.iter()
            .map(|p| p.distance_to(&pts[0]))
            .fold(0.0, f64::max)
    };
    let lin = tolerance.linear();
    let u_closed = (0..=8).all(|k| {
        let v = v0 + (v1 - v0) * k as f64 / 8.0;
        surface
            .point_at(u0, v)
            .distance_to(&surface.point_at(u1, v))
            <= lin
    });
    if u_closed {
        let pole = |v: f64| iso_spread(None, Some(v)) <= lin;
        if pole(v0) && pole(v1) {
            return replay_primitive_nurbs_ball::ball(shape, &surface, tolerance);
        }
        return closed_sheet(shape, &surface, tolerance);
    }
    for (name, spread) in [
        ("u start", iso_spread(Some(u0), None)),
        ("u end", iso_spread(Some(u1), None)),
        ("v start", iso_spread(None, Some(v0))),
        ("v end", iso_spread(None, Some(v1))),
    ] {
        if spread <= lin {
            return Err(invalid(
                shape,
                "points",
                spread,
                &format!(
                    "the net's {name} edge collapses to a point; only open patches are sheets"
                ),
            ));
        }
    }
    let mut brep = popkern_sweep::topology::side_faces::build_sweep_topology_with(
        std::slice::from_ref(&surface),
        &[],
        &[],
        false,
        popkern_sweep::topology::side_faces::SideFaceGeometry::Nurbs,
        tolerance,
        None,
    )
    .map_err(|e| failed(shape, "sheet", e))?;
    // The builder backs each free edge with a surfaceless back face;
    // the canonical sheet keeps the front alone, its boundary edges
    // lone fins, as every sheet primitive does.
    let back: Vec<_> = brep
        .faces
        .iter()
        .filter(|(_, f)| f.surface.is_none())
        .map(|(h, _)| h)
        .collect();
    if !back.is_empty() {
        brep.retire_faces_to_boundary(&back)
            .map_err(|e| failed(shape, "sheet", e))?;
    }
    brep.settle_open_shells();
    if brep.faces.len() != 1 {
        return Err(failed(
            shape,
            "sheet",
            format!(
                "the surface came as {} faces; one was expected",
                brep.faces.len()
            ),
        ));
    }
    adopt_mirror_loops(&mut brep);
    attach_boundary_curves(&mut brep, &surface, tolerance);
    Ok(brep)
}

/// The surface the net, weights and knots define, checked.
fn surface(
    shape: &str,
    points: &[Vec<[f64; 3]>],
    weights: &[Vec<f64>],
    degree: [u32; 2],
    knots: &[Vec<f64>; 2],
) -> Result<NurbsSurface, ReplayError> {
    let num_u = points.len();
    let num_v = points.first().map_or(0, Vec::len);
    let (du, dv) = (degree[0] as usize, degree[1] as usize);
    if du == 0 || dv == 0 {
        return Err(invalid(shape, "degree", 0.0, "degrees must be at least 1"));
    }
    if num_u < du + 1 || num_v < dv + 1 || points.iter().any(|r| r.len() != num_v) {
        return Err(invalid(
            shape,
            "points",
            num_u as f64,
            "needs a rectangular net of at least degree + 1 points each way",
        ));
    }
    let flat_weights: Vec<f64> = if weights.is_empty() {
        vec![1.0; num_u * num_v]
    } else {
        if weights.len() != num_u || weights.iter().any(|r| r.len() != num_v) {
            return Err(invalid(
                shape,
                "weights",
                weights.len() as f64,
                "must match the net",
            ));
        }
        weights.iter().flatten().copied().collect()
    };
    if let Some(w) = flat_weights
        .iter()
        .find(|w| w.is_nan() || **w <= 0.0 || !w.is_finite())
    {
        return Err(invalid(shape, "weights", *w, "must be positive"));
    }
    let knot_vector =
        |given: &Vec<f64>, n: usize, d: usize, name: &str| -> Result<Vec<f64>, ReplayError> {
            if given.is_empty() {
                // Clamped uniform.
                let inner = n - d - 1;
                let mut k = vec![0.0; d + 1];
                k.extend((1..=inner).map(|i| i as f64 / (inner + 1) as f64));
                k.extend(vec![1.0; d + 1]);
                return Ok(k);
            }
            if given.len() != n + d + 1 {
                return Err(invalid(
                    shape,
                    name,
                    given.len() as f64,
                    "needs points + degree + 1 knots",
                ));
            }
            if given.windows(2).any(|w| w[1] < w[0]) || given[0] >= given[given.len() - 1] {
                return Err(invalid(
                    shape,
                    name,
                    given[0],
                    "must be non-decreasing over a span",
                ));
            }
            let clamped = given[..=d].iter().all(|k| *k == given[0])
                && given[given.len() - d - 1..]
                    .iter()
                    .all(|k| *k == given[given.len() - 1]);
            if !clamped {
                return Err(invalid(
                    shape,
                    name,
                    given[0],
                    "must be clamped (end knots of multiplicity degree + 1)",
                ));
            }
            Ok(given.clone())
        };
    let knots_u = knot_vector(&knots[0], num_u, du, "knots_u")?;
    let knots_v = knot_vector(&knots[1], num_v, dv, "knots_v")?;
    Ok(NurbsSurface {
        degree_u: du,
        degree_v: dv,
        knots_u,
        knots_v,
        control_points: points
            .iter()
            .flatten()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect(),
        num_u,
        num_v,
        weights: flat_weights,
    })
}

/// The free boundary as lone fins, the canonical sheet edge (PK-978):
/// the builder leaves each boundary edge's other half-edge dangling, in
/// no loop, which no sheet consumer reads as open (the partition took
/// the sheet for a solid of unknown material sense). Each dangling twin
/// goes, re-anchored, and its edge keeps the face's half-edge alone.
fn adopt_mirror_loops(brep: &mut BRep) {
    use slotmap::Key;
    let dangling: Vec<_> = brep
        .half_edges
        .iter()
        .filter(|(_, h)| h.loop_.is_null() || !brep.loops.contains_key(h.loop_))
        .map(|(handle, _)| handle)
        .collect();
    for twin in dangling {
        let Some(survivor) = brep.half_edges[twin].mate() else {
            continue;
        };
        let edge = brep.half_edges[twin].edge;
        brep.remove_half_edge_reanchored(twin);
        if let Some(h) = brep.half_edges.get_mut(survivor) {
            h.clear_mate();
        }
        if let Some(e) = brep.edges.get_mut(edge) {
            e.half_edge = survivor;
        }
    }
}

/// The four boundary iso-curves of `surface`, each its net's boundary
/// row or column with their weights and the other direction's knots.
pub(super) fn boundary_curves(surface: &NurbsSurface) -> [NurbsCurve; 4] {
    let (nu, nv) = (surface.num_u, surface.num_v);
    let at = |i: usize, j: usize| i * nv + j;
    let row = |i: usize| NurbsCurve {
        degree: surface.degree_v,
        knots: surface.knots_v.clone(),
        control_points: (0..nv).map(|j| surface.control_points[at(i, j)]).collect(),
        weights: (0..nv).map(|j| surface.weights[at(i, j)]).collect(),
    };
    let column = |j: usize| NurbsCurve {
        degree: surface.degree_u,
        knots: surface.knots_u.clone(),
        control_points: (0..nu).map(|i| surface.control_points[at(i, j)]).collect(),
        weights: (0..nu).map(|i| surface.weights[at(i, j)]).collect(),
    };
    [row(0), row(nu - 1), column(0), column(nv - 1)]
}

/// Give each boundary edge its exact curve (the builder leaves an open
/// profile's boundary edges without one, and every consumer that walks
/// an edge, the edge-in-face recovery among them, skipped them): the
/// boundary iso-curve whose ends are the edge's vertices, run from the
/// edge's half-edge origin.
fn attach_boundary_curves(brep: &mut BRep, surface: &NurbsSurface, tolerance: Tolerance) {
    use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
    let curves = boundary_curves(surface);
    let reach = 10.0 * tolerance.linear();
    let edges: Vec<_> = brep.edges.keys().collect();
    for eh in edges {
        if brep.edges[eh].curve.is_some() {
            continue;
        }
        let he = brep.edges[eh].half_edge;
        let Some(h) = brep.half_edges.get(he) else {
            continue;
        };
        let far = h.mate().unwrap_or(h.next);
        let (Some(p), Some(q)) = (
            brep.vertices.get(h.origin).map(|v| v.position),
            brep.half_edges
                .get(far)
                .and_then(|f| brep.vertices.get(f.origin))
                .map(|v| v.position),
        ) else {
            continue;
        };
        let (u0, u1) = (
            surface.knots_u[0],
            surface.knots_u[surface.knots_u.len() - 1],
        );
        let (v0, v1) = (
            surface.knots_v[0],
            surface.knots_v[surface.knots_v.len() - 1],
        );
        for (k, c) in curves.iter().enumerate() {
            let (t0, t1) = c.parameter_range();
            let (a, b) = (c.point_at(t0), c.point_at(t1));
            let (t_start, t_end) = if a.distance_to(&p) <= reach && b.distance_to(&q) <= reach {
                (t0, t1)
            } else if b.distance_to(&p) <= reach && a.distance_to(&q) <= reach {
                (t1, t0)
            } else {
                continue;
            };
            brep.edges[eh].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(c.clone()),
                t_start,
                t_end,
            });
            // Its pcurve on the face: the iso-line it is in the chart.
            let uv = |t: f64| match k {
                0 => (u0, t),
                1 => (u1, t),
                2 => (t, v0),
                _ => (t, v1),
            };
            brep.half_edges[he].pcurve =
                Some(Pcurve::linear(uv(t_start), uv(t_end), t_start, t_end));
            break;
        }
    }
}

/// A surface closed round u (its first and last rows one curve): the
/// solid it bounds with flat discs at its rims (a B-spline cylinder's
/// side and its ends; a rim collapsed to a point, a cone's apex, takes a
/// degenerate cap), its seam the u = 0 iso-curve.
fn closed_sheet(
    shape: &str,
    surface: &NurbsSurface,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let rows = boundary_curves(surface);
    // rows: u start, u end, v start, v end (the last two run along u)
    let brep = popkern_sweep::topology::side_faces::build_sweep_topology_with(
        std::slice::from_ref(surface),
        std::slice::from_ref(&rows[2]),
        std::slice::from_ref(&rows[3]),
        true,
        popkern_sweep::topology::side_faces::SideFaceGeometry::Nurbs,
        tolerance,
        None,
    )
    .map_err(|e| failed(shape, "sheet", e))?;
    Ok(brep)
}
