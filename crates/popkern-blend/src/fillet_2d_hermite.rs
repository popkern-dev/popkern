//! Hermite corner blends for sheet-body wires (PK-508): a corner of a
//! planar wire replaced by a polynomial curve matching each support's
//! jets through order `k` — G1, G2 or G3 — instead of a G1 arc.
//!
//! The blend is the degree `n = 2k + 1` Bezier curve whose end control
//! points reproduce the supports' derivatives at the contact points:
//! `P1 = P0 + j1/n`, `P2 = 2 P1 - P0 + j2/(n (n-1))`, `P3 = 3 P2 - 3 P1
//! + P0 + j3/(n (n-1) (n-2))`, mirrored at the far end. Its parameters
//! are SETBACKS and SPEEDS, not a radius: each contact lies `setback`
//! along its support (arc length) from the corner, and the blend leaves
//! it at `speed` (its derivative's length there, per unit of its own
//! `[0, 1]` parameter); the supports' higher derivatives are those of
//! the support reparameterized to that speed. Rational supports (arcs,
//! NURBS) give their exact jets (`curve_join::JoinJet`).
//!
//! Before any surgery the blend is checked over its WHOLE length:
//! regular (its hodograph's control points all lie on one side of a
//! line through the origin, so its derivative never vanishes), free of
//! self-intersection, and clear of both supports away from its contacts.
//! A blend that fails declines naming the check.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::corner_2d_helpers::topology_failure;
use crate::curve_join::JoinJet;
use crate::error::BlendError;

/// A Hermite corner blend's parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HermiteCorner2d {
    /// Continuity order matched at both contacts: 1, 2 or 3.
    pub order: u8,
    /// Arc length from the corner to the contact on the first support
    /// (the support whose outgoing half-edge sorts first).
    pub setback_a: f64,
    /// The same on the second support.
    pub setback_b: f64,
    /// The blend's speed leaving the first contact.
    pub speed_a: f64,
    /// Its speed arriving at the second.
    pub speed_b: f64,
}

/// Replace the corner at `vertex` of a sheet body's wire by the Hermite
/// blend `blend` (see the module docs).
///
/// # Errors
///
/// [`BlendError::InvalidParameter`] for an order outside 1..=3 or a
/// non-positive setback or speed; [`BlendError::InvalidTopology`] when
/// the corner is not a two-edge wire corner, a setback reaches past its
/// support, or the blend is irregular, self-intersecting or crosses a
/// support.
pub fn hermite_corner_2d(
    brep: &BRep,
    vertex: VertexHandle,
    blend: &HermiteCorner2d,
) -> Result<BRep, BlendError> {
    let bad = |parameter: &str, value: f64, reason: &str| BlendError::InvalidParameter {
        parameter: parameter.into(),
        value: format!("{value}"),
        reason: reason.into(),
    };
    if !(1..=3).contains(&blend.order) {
        return Err(bad(
            "order",
            f64::from(blend.order),
            "1 (G1), 2 (G2) or 3 (G3)",
        ));
    }
    for (name, v) in [
        ("setback_a", blend.setback_a),
        ("setback_b", blend.setback_b),
        ("speed_a", blend.speed_a),
        ("speed_b", blend.speed_b),
    ] {
        if !(v.is_finite() && v > 0.0) {
            return Err(bad(name, v, "a positive finite length"));
        }
    }
    let curve =
        plan(brep, vertex, blend).map_err(|reason| topology_failure(brep, vertex, reason))?;
    let (he_a, he_b) =
        outgoing_pair(brep, vertex).map_err(|r| topology_failure(brep, vertex, r))?;
    let plan = crate::fillet_2d_plan::CornerPlan::for_curve(he_a, he_b, curve);
    let mut out = brep.clone();
    crate::fillet_2d::apply_corner_plan(&mut out, vertex, &plan)
        .map_err(|reason| topology_failure(brep, vertex, reason))?;
    out.attach_planar_pcurves()
        .map_err(|e| topology_failure(brep, vertex, format!("pcurve refresh failed: {e:?}")))?;
    out.validate()
        .map_err(|e| topology_failure(brep, vertex, format!("validation failed: {e:?}")))?;
    Ok(out)
}

/// The corner's two outgoing half-edges, sorted (as `corner_plan`).
fn outgoing_pair(
    brep: &BRep,
    vertex: VertexHandle,
) -> Result<(HalfEdgeHandle, HalfEdgeHandle), String> {
    let mut out: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == vertex)
        .map(|(h, _)| h)
        .collect();
    out.sort();
    match out.as_slice() {
        [a, b] if brep.half_edges[*a].edge != brep.half_edges[*b].edge => Ok((*a, *b)),
        _ => Err(format!(
            "corner vertex must join two edges by two outgoing half-edges, got {}",
            out.len()
        )),
    }
}

/// A support of the corner: its curve, the parameter at the contact,
/// and the sign of the parameter direction running TOWARD the corner.
struct Contact {
    jet: JoinJet,
    toward_corner: f64,
}

fn contact(
    brep: &BRep,
    he: HalfEdgeHandle,
    corner: Point3,
    setback: f64,
) -> Result<Contact, String> {
    let edge = &brep.edges[brep.half_edges[he].edge];
    let curve = edge.curve.as_ref().ok_or("a support edge has no curve")?;
    let (t0, t1) = (curve.t_start, curve.t_end);
    let at_start = curve.geometry.point_at(t0).distance_to(&corner)
        <= curve.geometry.point_at(t1).distance_to(&corner);
    let total = curve.geometry.length_between(t0.min(t1), t0.max(t1));
    if setback >= total {
        return Err(format!(
            "setback {setback:.6} reaches past its support (length {total:.6})"
        ));
    }
    // The parameter `setback` along the support from the corner.
    let from = if at_start { t0 } else { t1 };
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let t = from + (if at_start { t1 - t0 } else { t0 - t1 }) * mid;
        let s = curve.geometry.length_between(from.min(t), from.max(t));
        if s < setback {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t = from + (if at_start { t1 - t0 } else { t0 - t1 }) * 0.5 * (lo + hi);
    let toward = if at_start == (t1 > t0) { -1.0 } else { 1.0 };
    Ok(Contact {
        jet: JoinJet::of(&curve.geometry, t),
        toward_corner: toward,
    })
}

/// The blend curve, checked over its whole length.
fn plan(brep: &BRep, vertex: VertexHandle, blend: &HermiteCorner2d) -> Result<NurbsCurve, String> {
    let (he_a, he_b) = outgoing_pair(brep, vertex)?;
    let corner = brep.vertices[vertex].position;
    let a = contact(brep, he_a, corner, blend.setback_a)?;
    let b = contact(brep, he_b, corner, blend.setback_b)?;
    let k = usize::from(blend.order);
    let n = 2 * k + 1;
    // Jets of each support reparameterized to the blend's speed: on `a`
    // running toward the corner, on `b` running away from it.
    let scaled = |c: &Contact, speed: f64, sense: f64| -> Result<[Vec3; 3], String> {
        let len = c.jet.d[0].length();
        if len <= 0.0 {
            return Err("a support is singular at its contact".into());
        }
        let lam = sense * speed / len;
        Ok([
            c.jet.d[0] * lam,
            c.jet.d[1] * (lam * lam),
            c.jet.d[2] * lam.powi(3),
        ])
    };
    let j = scaled(&a, blend.speed_a, a.toward_corner)?;
    let e = scaled(&b, blend.speed_b, -b.toward_corner)?;
    let nf = n as f64;
    let mut p = vec![Point3::new(0.0, 0.0, 0.0); n + 1];
    p[0] = a.jet.point;
    p[n] = b.jet.point;
    // Start: P1 = P0 + j1/n, P2 = 2P1 - P0 + j2/(n(n-1)), P3 = ...
    if k >= 1 {
        p[1] = p[0] + j[0] * (1.0 / nf);
        p[n - 1] = p[n] - e[0] * (1.0 / nf);
    }
    if k >= 2 {
        let f = 1.0 / (nf * (nf - 1.0));
        p[2] = p[1] + (p[1] - p[0]) + j[1] * f;
        p[n - 2] = p[n - 1] + (p[n - 1] - p[n]) + e[1] * f;
    }
    if k >= 3 {
        let f = 1.0 / (nf * (nf - 1.0) * (nf - 2.0));
        let v = |q: Point3| Vec3::new(q.x, q.y, q.z);
        let d = v(p[2]) * 3.0 - v(p[1]) * 3.0 + v(p[0]) + j[2] * f;
        p[3] = Point3::new(d.x, d.y, d.z);
        let d = v(p[n - 2]) * 3.0 - v(p[n - 1]) * 3.0 + v(p[n]) - e[2] * f;
        p[n - 3] = Point3::new(d.x, d.y, d.z);
    }
    let mut knots = vec![0.0; n + 1];
    knots.extend(vec![1.0; n + 1]);
    let curve = NurbsCurve {
        degree: n,
        knots,
        control_points: p.clone(),
        weights: vec![1.0; n + 1],
    };
    check_regular(&p)?;
    let samples: Vec<Point3> = (0..=256)
        .map(|i| curve.point_at(i as f64 / 256.0))
        .collect();
    check_simple(&samples)?;
    for (he, end) in [(he_a, 0usize), (he_b, samples.len() - 1)] {
        check_clear(brep, he, &samples, end)?;
    }
    Ok(curve)
}

/// The hodograph's control points `n (P_{i+1} - P_i)` all on one side
/// of a line through the origin: the derivative's convex hull misses
/// zero, so it never vanishes.
fn check_regular(p: &[Point3]) -> Result<(), String> {
    let h: Vec<Vec3> = p.windows(2).map(|w| w[1] - w[0]).collect();
    let mean = h.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, v| s + *v);
    let dir = mean
        .normalized()
        .ok_or("the blend's hodograph sums to zero")?;
    if h.iter().all(|v| v.dot(&dir) > 0.0) {
        Ok(())
    } else {
        Err(
            "the blend is not certified regular (its hodograph hull reaches the origin): \
             lower the speeds or move the setbacks"
                .into(),
        )
    }
}

fn segments_cross(a: Point3, b: Point3, c: Point3, d: Point3) -> bool {
    // The wire is planar: test in the plane of the four points' first
    // two independent directions via 3-D orientation signs.
    let n = (b - a).cross(&(d - c));
    let n = if n.length() > 1e-300 {
        n
    } else {
        return false;
    };
    let side = |o: Point3, p: Point3, q: Point3| (p - o).cross(&(q - o)).dot(&n);
    side(a, b, c) * side(a, b, d) < 0.0 && side(c, d, a) * side(c, d, b) < 0.0
}

fn check_simple(s: &[Point3]) -> Result<(), String> {
    for i in 0..s.len() - 1 {
        for j in i + 2..s.len() - 1 {
            if segments_cross(s[i], s[i + 1], s[j], s[j + 1]) {
                return Err("the blend crosses itself".into());
            }
        }
    }
    Ok(())
}

/// The blend must not cross the kept part of the support of half-edge
/// `he` (sampled), away from its own contact sample `end`.
fn check_clear(
    brep: &BRep,
    he: HalfEdgeHandle,
    blend: &[Point3],
    end: usize,
) -> Result<(), String> {
    let edge = &brep.edges[brep.half_edges[he].edge];
    let Some(curve) = edge.curve.as_ref() else {
        return Ok(());
    };
    let support: Vec<Point3> = (0..=128)
        .map(|i| {
            curve
                .geometry
                .point_at(curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 128.0)
        })
        .collect();
    for (i, w) in blend.windows(2).enumerate() {
        // The segments touching the contact meet the support there.
        if i + 1 == end || i == end || i + 2 == end || i == end + 1 {
            continue;
        }
        for v in support.windows(2) {
            if segments_cross(w[0], w[1], v[0], v[1]) {
                return Err("the blend crosses a support".into());
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "fillet_2d_hermite_tests.rs"]
mod tests;
