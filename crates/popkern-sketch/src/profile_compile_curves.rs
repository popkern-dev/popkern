//! Conic and spline pieces of an analytic profile (split from
//! `profile_compile` under the file budget): validation, their
//! entries in the geometry-container sketch, and the general wire
//! assembly a path with a spline segment needs.
//!
//! Everything stays exact: a conic is one rational quadratic span, an
//! ellipse four, a spline its own NURBS. A path of lines, arcs and
//! conics keeps the quadratic span assembly; a path with a spline is
//! every piece raised to the highest degree and joined C0, each piece
//! on its own unit knot interval.

use popkern_constraint::{PointHandle, Sketch};
use popkern_core::geometry::{NurbsCurve, Point3};

use super::{Piece, arc_run};
use crate::exact_arc::{QuadraticSpan, assemble_quadratic, line_span};
use crate::profile::ProfileError;

/// How far apart two pieces' shared ends may be when joined: they are
/// the same description coordinates, so only round-off separates them.
const JOIN_GAP: f64 = 1.0e-9;

fn invalid(what: &'static str, reason: impl Into<String>) -> ProfileError {
    ProfileError::InvalidSpline {
        what,
        reason: reason.into(),
    }
}

fn point(p: (f64, f64)) -> Point3 {
    Point3::new(p.0, p.1, 0.0)
}

/// Checked weights for `n` control points (empty: all 1).
fn weights(what: &'static str, weights: &[f64], n: usize) -> Result<Vec<f64>, ProfileError> {
    if weights.is_empty() {
        return Ok(vec![1.0; n]);
    }
    if weights.len() != n {
        return Err(invalid(
            what,
            format!("{} weights for {n} control points", weights.len()),
        ));
    }
    if let Some(w) = weights.iter().find(|w| !(w.is_finite() && **w > 0.0)) {
        return Err(invalid(what, format!("weight {w} is not positive")));
    }
    Ok(weights.to_vec())
}

fn degree(what: &'static str, degree: usize, n: usize) -> Result<usize, ProfileError> {
    if degree == 0 || degree >= n {
        return Err(invalid(
            what,
            format!(
                "degree {degree} needs between 1 and {} for {n} control points",
                n.max(2) - 1
            ),
        ));
    }
    Ok(degree)
}

/// The conic from `a` to `b` with shoulder `control` and middle
/// weight `w`.
pub(super) fn conic_piece(
    a: (f64, f64),
    control: (f64, f64),
    b: (f64, f64),
    w: f64,
) -> Result<Piece, ProfileError> {
    if !(w.is_finite() && w > 0.0) {
        return Err(invalid("conic_to", format!("weight {w} is not positive")));
    }
    Ok(Piece::Quadratic(QuadraticSpan {
        p0: point(a),
        pm: point(control),
        p2: point(b),
        w,
        t0: 0.0,
        t1: 1.0,
    }))
}

/// The four exact quarters of the ellipse of semi-axes `rx`, `ry`
/// about `(cx, cy)`, counter-clockwise from `+x`.
pub(super) fn ellipse_pieces(cx: f64, cy: f64, rx: f64, ry: f64) -> Vec<Piece> {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let on = |x: f64, y: f64| point((cx + x, cy + y));
    let ends = [(rx, 0.0), (0.0, ry), (-rx, 0.0), (0.0, -ry)];
    (0..4)
        .map(|i| {
            let (a, b) = (ends[i], ends[(i + 1) % 4]);
            Piece::Quadratic(QuadraticSpan {
                p0: on(a.0, a.1),
                // The shoulder: the corner of the bounding box the two
                // ends' tangents meet at.
                pm: on(a.0 + b.0, a.1 + b.1),
                p2: on(b.0, b.1),
                w,
                t0: i as f64 / 4.0,
                t1: (i + 1) as f64 / 4.0,
            })
        })
        .collect()
}

/// The clamped NURBS from `a` through `control` to `b`.
pub(super) fn spline_piece(
    a: (f64, f64),
    control: &[(f64, f64)],
    b: (f64, f64),
    p: usize,
    ws: &[f64],
    knots: &[f64],
) -> Result<Piece, ProfileError> {
    const WHAT: &str = "spline_to";
    let mut control_points = vec![point(a)];
    control_points.extend(control.iter().copied().map(point));
    control_points.push(point(b));
    let n = control_points.len();
    let p = degree(WHAT, p, n)?;
    let weights = weights(WHAT, ws, n)?;
    let knots = if knots.is_empty() {
        clamped_uniform(n, p)
    } else {
        clamped_knots(WHAT, knots, n, p)?
    };
    Ok(Piece::Nurbs(NurbsCurve {
        degree: p,
        knots,
        control_points,
        weights,
    }))
}

/// The closed periodic B-spline of control polygon `points` (offset
/// already applied), uniform knots.
pub(super) fn periodic_piece(
    points: Vec<(f64, f64)>,
    p: usize,
    ws: &[f64],
) -> Result<Piece, ProfileError> {
    const WHAT: &str = "spline";
    let n = points.len();
    if n < 3 {
        return Err(invalid(
            WHAT,
            format!("a closed spline needs at least 3 control points, got {n}"),
        ));
    }
    let degree = degree(WHAT, p, n)?;
    let weights = weights(WHAT, ws, n)?;
    Ok(Piece::Periodic {
        points,
        degree,
        weights,
    })
}

fn clamped_uniform(n: usize, p: usize) -> Vec<f64> {
    let spans = (n - p) as f64;
    let mut k = vec![0.0; p + 1];
    k.extend((1..n - p).map(|i| i as f64 / spans));
    k.extend(std::iter::repeat_n(1.0, p + 1));
    k
}

/// A user's knot vector, checked: `n + p + 1` long, non-decreasing,
/// clamped at both ends, no interior knot past multiplicity `p`.
fn clamped_knots(
    what: &'static str,
    knots: &[f64],
    n: usize,
    p: usize,
) -> Result<Vec<f64>, ProfileError> {
    if knots.len() != n + p + 1 {
        return Err(invalid(
            what,
            format!(
                "{} knots for {n} control points of degree {p} (needs {})",
                knots.len(),
                n + p + 1
            ),
        ));
    }
    if knots.iter().any(|k| !k.is_finite()) || knots.windows(2).any(|w| w[1] < w[0]) {
        return Err(invalid(what, "knots must be finite and non-decreasing"));
    }
    let (first, last) = (knots[0], knots[knots.len() - 1]);
    if last <= first
        || knots[..=p].iter().any(|k| *k != first)
        || knots[knots.len() - p - 1..].iter().any(|k| *k != last)
    {
        return Err(invalid(
            what,
            format!("knots must be clamped: the first and last {} equal", p + 1),
        ));
    }
    let interior = &knots[p + 1..knots.len() - p - 1];
    for k in interior {
        if interior.iter().filter(|x| *x == k).count() > p {
            return Err(invalid(
                what,
                format!("interior knot {k} repeats past degree {p}"),
            ));
        }
    }
    Ok(knots.to_vec())
}

/// The periodic spline as the one-period unclamped NURBS the
/// constraint dialect evaluates (extended knots, wrapped control
/// points), parameter range `[0, n]`.
fn periodic_unclamped(points: &[(f64, f64)], p: usize, weights: &[f64]) -> NurbsCurve {
    let n = points.len();
    NurbsCurve {
        degree: p,
        knots: (0..n + 2 * p + 1).map(|i| i as f64 - p as f64).collect(),
        control_points: (0..n + p).map(|j| point(points[j % n])).collect(),
        weights: (0..n + p).map(|j| weights[j % n]).collect(),
    }
}

/// `curve` over `[a, b]` in clamped form: each end knot inserted to
/// multiplicity `p`, so the curve passes through a control point
/// there, then the window cut out with its ends clamped. The same
/// curve, exactly (knot insertion does not move it).
fn clamp_window(curve: &NurbsCurve, a: f64, b: f64) -> Result<NurbsCurve, ProfileError> {
    let p = curve.degree;
    let mut c = curve.clone();
    for t in [a, b] {
        let have = c.knots.iter().filter(|k| **k == t).count();
        for _ in have..p {
            c = c
                .insert_knot(t)
                .map_err(|e| invalid("spline", format!("clamping the closed spline: {e}")))?;
        }
    }
    let first = |t: f64| c.knots.iter().position(|k| *k == t);
    let (Some(s), Some(e)) = (first(a), first(b)) else {
        return Err(invalid(
            "spline",
            "clamping the closed spline lost its window",
        ));
    };
    let mut knots = vec![a; p + 1];
    knots.extend_from_slice(&c.knots[s + p..e]);
    knots.extend(std::iter::repeat_n(b, p + 1));
    Ok(NurbsCurve {
        degree: p,
        knots,
        control_points: c.control_points[s - 1..e].to_vec(),
        weights: c.weights[s - 1..e].to_vec(),
    })
}

/// The closed periodic spline as one clamped closed wire.
pub(super) fn periodic_wire(
    points: &[(f64, f64)],
    p: usize,
    weights: &[f64],
) -> Result<NurbsCurve, ProfileError> {
    clamp_window(
        &periodic_unclamped(points, p, weights),
        0.0,
        points.len() as f64,
    )
}

/// One piece as a clamped NURBS on the knot interval `[0, 1]`.
fn piece_curve(piece: &Piece) -> Result<NurbsCurve, ProfileError> {
    let mut c = match piece {
        Piece::Line(a, b) => assemble_quadratic(&[line_span(point(*a), point(*b), 0.0, 1.0)]),
        Piece::Arc(arc) => {
            // Uniform labels: a clockwise run's angle labels descend.
            let mut run = arc_run(arc);
            let n = run.len() as f64;
            for (i, s) in run.iter_mut().enumerate() {
                (s.t0, s.t1) = (i as f64 / n, (i + 1) as f64 / n);
            }
            assemble_quadratic(&run)
        }
        Piece::Quadratic(span) => assemble_quadratic(&[*span]),
        Piece::Nurbs(curve) => curve.clone(),
        // A closed spline is a whole shape, never a path piece.
        Piece::Periodic { .. } => {
            return Err(invalid("spline", "a closed spline is not a path segment"));
        }
    };
    let (k0, k1) = (c.knots[0], c.knots[c.knots.len() - 1]);
    for k in &mut c.knots {
        *k = (*k - k0) / (k1 - k0);
    }
    Ok(c)
}

/// A path's pieces joined into one closed wire: every piece raised to
/// the highest degree and concatenated C0.
pub(super) fn general_wire(pieces: &[Piece]) -> Result<NurbsCurve, ProfileError> {
    let curves = pieces
        .iter()
        .map(piece_curve)
        .collect::<Result<Vec<_>, _>>()?;
    let p = curves.iter().map(|c| c.degree).max().unwrap_or(1);
    let mut out: Option<NurbsCurve> = None;
    for c in curves {
        let c = if c.degree < p {
            c.elevate_degree_to(p)
                .map_err(|e| invalid("spline_to", format!("raising a piece's degree: {e}")))?
        } else {
            c
        };
        out = Some(match out {
            None => c,
            Some(acc) => acc
                .concatenate_clamped(&c, JOIN_GAP)
                .map_err(|e| invalid("spline_to", format!("joining the path: {e}")))?,
        });
    }
    out.ok_or(ProfileError::PathTooShort)
}

/// Load a conic or spline piece into the geometry-container sketch
/// (its own control points; the edge graph joins coincident ends).
pub(super) fn add_curve_element(sk: &mut Sketch, piece: &Piece) {
    let mut add = |cps: &[Point3]| -> Vec<PointHandle> {
        cps.iter().map(|c| sk.add_point(c.x, c.y)).collect()
    };
    match piece {
        Piece::Quadratic(span) => {
            let h = add(&[span.p0, span.pm, span.p2]);
            let s = sk.add_spline(&h, 2);
            let _ = sk.set_spline_weight(s, 1, span.w);
        }
        Piece::Nurbs(curve) => {
            let h = add(&curve.control_points);
            let s = sk.add_spline(&h, curve.degree);
            for (i, w) in curve.weights.iter().enumerate() {
                let _ = sk.set_spline_weight(s, i, *w);
            }
            for (i, k) in curve.knots.iter().enumerate() {
                let _ = sk.set_spline_knot(s, i, *k);
            }
        }
        Piece::Periodic {
            points,
            degree,
            weights,
        } => {
            let h: Vec<PointHandle> = points.iter().map(|q| sk.add_point(q.0, q.1)).collect();
            sk.add_periodic_spline(&h, *degree, weights);
        }
        Piece::Line(..) | Piece::Arc(_) => {}
    }
}

#[cfg(test)]
#[path = "profile_compile_curves_tests.rs"]
mod tests;
