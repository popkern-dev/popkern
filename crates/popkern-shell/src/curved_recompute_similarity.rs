//! Similarity-map recompute arms for curved inner-skin edges (PK-324).
//!
//! Two arms, both expressing the recomputed curve as an exact
//! similarity image `q = A + s·p` of the clone curve — the form every
//! analytic offset of a line/circle edge takes — so parametrization and
//! direction carry over by construction (the PK-303 direction contract):
//!
//! - **Parallel-map arm**: seam edges and G1-tangential pairs. Along
//!   such an edge both offset surfaces displace each point identically
//!   (they share the surface normal there), so the exact closest-point
//!   projection onto either offset surface IS the normal-offset map.
//!   The agreement of the two projections is the class witness — no
//!   pair enumeration, measured membership.
//! - **Coaxial-profile arm**: transversal pairs of revolution surfaces
//!   about a common axis (bore×countersink cylinder×cone junctions).
//!   Intersection circles are axis-centered stations where the two
//!   profile radius functions agree; solved in the 1D profile plane,
//!   Newton seeded from the clone circle's own station.
//!
//! Every candidate must pass through the recomputed vertices and lie on
//! both offset surfaces at sampled parameters, else `None` — the clone
//! stays stale and the serve gate declines the body.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;

use crate::curved_recompute::{analytic_residual, project_to_analytic_surface};

#[path = "curved_recompute_profile.rs"]
mod curved_recompute_profile;
pub(crate) use curved_recompute_profile::coaxial_profile_edge_curve;

/// Map a curve geometry through the similarity `p ↦ anchor + s·(p−0)`,
/// preserving parametrization and direction exactly. Circle/Arc frames
/// derive from the (unchanged) normal, so evaluation commutes with the
/// map; rational NURBS evaluation commutes with any affine map of the
/// control points.
pub(crate) fn map_similarity_geometry(
    geometry: &CurveGeometry,
    anchor: Vec3,
    s: f64,
) -> Option<CurveGeometry> {
    let map_pt =
        |p: &Point3| Point3::new(anchor.x + s * p.x, anchor.y + s * p.y, anchor.z + s * p.z);
    match geometry {
        CurveGeometry::Line { origin, direction } => Some(CurveGeometry::Line {
            origin: map_pt(origin),
            direction: *direction * s,
        }),
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => Some(CurveGeometry::Circle {
            center: map_pt(center),
            normal: *normal,
            radius: radius * s,
        }),
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        } => Some(CurveGeometry::Arc {
            center: map_pt(center),
            normal: *normal,
            radius: radius * s,
            start_angle: *start_angle,
            end_angle: *end_angle,
        }),
        CurveGeometry::Nurbs(n) => {
            let mut mapped = n.clone();
            for cp in &mut mapped.control_points {
                *cp = map_pt(cp);
            }
            Some(CurveGeometry::Nurbs(mapped))
        }
        _ => None,
    }
}

const SAMPLE_FRACTIONS: [f64; 7] = [0.0, 1.0 / 6.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 5.0 / 6.0, 1.0];

/// Verify a candidate: endpoints at the recomputed vertex positions and
/// sampled residuals on both offset surfaces within the linear
/// tolerance.
pub(super) fn verify_candidate(
    candidate: &CurveOnEdge,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    start_pos: Point3,
    end_pos: Point3,
    lin: f64,
) -> bool {
    if candidate
        .geometry
        .point_at(candidate.t_start)
        .distance_to(&start_pos)
        > lin
        || candidate
            .geometry
            .point_at(candidate.t_end)
            .distance_to(&end_pos)
            > lin
    {
        return false;
    }
    let (t0, t1) = (candidate.t_start, candidate.t_end);
    for f in SAMPLE_FRACTIONS {
        let p = candidate.geometry.point_at(t0 + (t1 - t0) * f);
        for surf in [surf1, surf2] {
            if analytic_residual(surf, &p).is_none_or(|r| r > lin) {
                return false;
            }
        }
    }
    true
}

/// Fit the similarity `q = A + s·p` from paired samples and map the
/// clone geometry through it. `None` when the samples are degenerate or
/// don't actually follow one similarity.
fn similarity_from_pairs(
    old: &CurveOnEdge,
    pairs: &[(Point3, Point3)],
    lin: f64,
) -> Option<CurveGeometry> {
    // Widest sample pair sets the scale.
    let (mut bi, mut bj, mut best) = (0usize, 0usize, 0.0_f64);
    for i in 0..pairs.len() {
        for j in i + 1..pairs.len() {
            let d = pairs[i].0.distance_to(&pairs[j].0);
            if d > best {
                (bi, bj, best) = (i, j, d);
            }
        }
    }
    if best <= lin {
        return None;
    }
    let s = pairs[bi].1.distance_to(&pairs[bj].1) / best;
    if !(s.is_finite() && s > 0.0) {
        return None;
    }
    let a = Vec3::new(
        pairs[bi].1.x - s * pairs[bi].0.x,
        pairs[bi].1.y - s * pairs[bi].0.y,
        pairs[bi].1.z - s * pairs[bi].0.z,
    );
    // Every projected sample must sit on the fitted similarity — a
    // rotation-carrying or non-uniform map is not this arm's class.
    for (p, q) in pairs {
        let mapped = Point3::new(a.x + s * p.x, a.y + s * p.y, a.z + s * p.z);
        if mapped.distance_to(q) > lin * 0.5 {
            return None;
        }
    }
    map_similarity_geometry(&old.geometry, a, s)
}

/// Parallel-map arm: seam edges (`surf1 == surf2`) and G1-tangential
/// pairs. Class witness: the exact projections of each clone sample
/// onto the two offset surfaces agree (they share the displacement
/// along the common normal). Transversal pairs fail the witness and
/// fall through to the other arms.
pub(crate) fn parallel_map_edge_curve(
    old: &CurveOnEdge,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    start_pos: Point3,
    end_pos: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    let lin = tolerance.linear();
    let (t0, t1) = (old.t_start, old.t_end);
    let trace = popkern_core::trace::trace_enabled("recompute");
    let mut pairs = Vec::with_capacity(SAMPLE_FRACTIONS.len());
    for f in SAMPLE_FRACTIONS {
        let p = old.geometry.point_at(t0 + (t1 - t0) * f);
        let Some(q1) = project_to_analytic_surface(surf1, &p) else {
            if trace {
                eprintln!("[recomp] pm: proj1 none at f={f}");
            }
            return None;
        };
        let Some(q2) = project_to_analytic_surface(surf2, &p) else {
            if trace {
                eprintln!("[recomp] pm: proj2 none at f={f}");
            }
            return None;
        };
        if q1.distance_to(&q2) > lin * 0.5 {
            if trace {
                eprintln!("[recomp] pm: disagree {} at f={f}", q1.distance_to(&q2));
            }
            return None; // transversal pair — not this arm's class
        }
        let q = Point3::new(
            0.5 * (q1.x + q2.x),
            0.5 * (q1.y + q2.y),
            0.5 * (q1.z + q2.z),
        );
        pairs.push((p, q));
    }
    let Some(geometry) = similarity_from_pairs(old, &pairs, lin) else {
        if trace {
            eprintln!("[recomp] pm: similarity fit declined");
        }
        return None;
    };
    let candidate = CurveOnEdge {
        geometry,
        t_start: t0,
        t_end: t1,
    };
    if verify_candidate(&candidate, surf1, surf2, start_pos, end_pos, lin) {
        return Some(candidate);
    }
    // Endpoints may have moved ALONG the curve: a seam arc ending on a
    // transversal junction keeps its geometry under the parallel map
    // but the junction vertex slides to a different station (measured:
    // sphere-capped rod, the sphere meridian's junction end moves
    // 3.5e-3 along the offset meridian). Re-trim to the recomputed
    // vertices, preserving traversal direction.
    if let Some(retrimmed) = retrim_to_endpoints(&candidate, start_pos, end_pos, lin)
        && verify_candidate(&retrimmed, surf1, surf2, start_pos, end_pos, lin)
    {
        return Some(retrimmed);
    }
    if trace {
        eprintln!(
            "[recomp] pm: verify declined (start {:?} cand {:?})",
            start_pos,
            candidate.geometry.point_at(candidate.t_start)
        );
    }
    None
}

/// Re-trim an open candidate so its parameter interval ends at the
/// recomputed vertex positions, preserving traversal direction. `None`
/// for closed edges (single vertex — the original interval is the
/// contract there) or when the vertices don't lie on the curve.
fn retrim_to_endpoints(
    candidate: &CurveOnEdge,
    start_pos: Point3,
    end_pos: Point3,
    lin: f64,
) -> Option<CurveOnEdge> {
    if start_pos.distance_to(&end_pos) <= lin {
        return None;
    }
    let ts = candidate.geometry.closest_point(&start_pos)?;
    let te = candidate.geometry.closest_point(&end_pos)?;
    if candidate.geometry.point_at(ts).distance_to(&start_pos) > lin
        || candidate.geometry.point_at(te).distance_to(&end_pos) > lin
    {
        return None;
    }
    if ((te - ts) > 0.0) != ((candidate.t_end - candidate.t_start) > 0.0) {
        return None;
    }
    Some(CurveOnEdge {
        geometry: candidate.geometry.clone(),
        t_start: ts,
        t_end: te,
    })
}

/// General numeric SSI arm (design §5 Stage 2): intersect the two
/// offset surfaces over windows bounded around the clone curve's own
/// samples, pick the branch nearest the clone midpoint, and re-trim to
/// the recomputed vertices. This is the arm for pairs with no analytic
/// form — NURBS offsets against anything. The candidate is verified
/// against the STORED offset surfaces (for approximated offsets that
/// is the serving contract: the body is consistent with its own
/// geometry; the approximation error is measured and reported by the
/// offset stage itself).
pub(crate) fn general_ssi_edge_curve(
    old: &CurveOnEdge,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    start_pos: Point3,
    end_pos: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    use popkern_core::geometry::intersection::surface_surface::{
        SurfaceSurfaceIntersection, intersect_surfaces,
    };
    use popkern_core::geometry::uv_projection::bounded_parameter_range;

    let lin = tolerance.linear();
    let (t0, t1) = (old.t_start, old.t_end);
    let samples: Vec<Point3> = SAMPLE_FRACTIONS
        .iter()
        .map(|f| old.geometry.point_at(t0 + (t1 - t0) * f))
        .collect();
    let (u1r, v1r) = bounded_parameter_range(&samples, surf1);
    let (u2r, v2r) = bounded_parameter_range(&samples, surf2);
    if ![u1r, v1r, u2r, v2r]
        .iter()
        .all(|(a, b)| a.is_finite() && b.is_finite())
    {
        return None;
    }

    let result = intersect_surfaces(surf1, u1r, v1r, surf2, u2r, v2r, tolerance).ok()?;
    // Branch selection: nearest to the clone curve's midpoint (the
    // center_curve pattern the design names).
    let mid = old.geometry.point_at(t0 + (t1 - t0) * 0.5);
    let mut best: Option<(f64, CurveGeometry)> = None;
    for intersection in &result.intersections {
        let SurfaceSurfaceIntersection::Curve(cos) = intersection else {
            continue;
        };
        let geometry = CurveGeometry::Nurbs(cos.curve_3d.clone());
        let Some(tm) = geometry.closest_point(&mid) else {
            continue;
        };
        let d = geometry.point_at(tm).distance_to(&mid);
        if best.as_ref().is_none_or(|(bd, _)| d < *bd) {
            best = Some((d, geometry));
        }
    }
    let (branch_distance, geometry) = best?;
    // The chosen branch must actually run near the clone curve — a
    // far branch means the true intersection near this edge was not
    // found and the clone must stay stale (structural decline).
    if branch_distance > 10.0 * lin + 0.5 * start_pos.distance_to(&end_pos) {
        return None;
    }

    // Trim to the recomputed vertices, direction from the vertex order.
    let ts = geometry.closest_point(&start_pos)?;
    let te = geometry.closest_point(&end_pos)?;
    if geometry.point_at(ts).distance_to(&start_pos) > lin
        || geometry.point_at(te).distance_to(&end_pos) > lin
        || (te - ts).abs() < 1.0e-12
    {
        return None;
    }
    let candidate = CurveOnEdge {
        geometry,
        t_start: ts,
        t_end: te,
    };
    // Residual verification against both stored offset surfaces.
    for f in SAMPLE_FRACTIONS {
        let p = candidate.geometry.point_at(ts + (te - ts) * f);
        for surf in [surf1, surf2] {
            if surface_residual(surf, &p).is_none_or(|r| r > lin) {
                return None;
            }
        }
    }
    Some(candidate)
}

/// Distance from a point to any surface: analytic closed forms where
/// they exist, iterative UV projection otherwise (NURBS offsets).
fn surface_residual(surface: &SurfaceGeometry, p: &Point3) -> Option<f64> {
    if let Some(r) = analytic_residual(surface, p) {
        return Some(r);
    }
    use popkern_core::geometry::surface::SurfaceEval;
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, p)?;
    Some(surface.point_at(u, v).distance_to(p))
}
