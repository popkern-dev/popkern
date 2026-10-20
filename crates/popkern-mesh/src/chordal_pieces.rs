//! Parameter pieces for the chordal certificate (see
//! `chordal_certificate`): seam unwrapping, pole fans, C⁰ breaks.

use popkern_core::geometry::Point3;
use popkern_core::geometry::certify::chord::C0Breaks;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Every vertex within half a period of one base vertex in each
/// periodic direction (a seam-straddling triangle otherwise reads as
/// spanning the whole period). The base giving the smallest span wins:
/// a pole vertex's parameter is arbitrary and must not be the base.
pub(crate) fn unwrap(w: [(f64, f64); 3], periods: (Option<f64>, Option<f64>)) -> [(f64, f64); 3] {
    let fold = |x: f64, base: f64, p: Option<f64>| match p.filter(|p| *p > 0.0) {
        Some(p) => x - ((x - base) / p).round() * p,
        None => x,
    };
    let around = |k: usize| -> [(f64, f64); 3] {
        let base = w[k];
        w.map(|q| (fold(q.0, base.0, periods.0), fold(q.1, base.1, periods.1)))
    };
    let span = |t: &[(f64, f64); 3]| {
        let ext = |f: fn(&(f64, f64)) -> f64| {
            let xs = t.iter().map(f);
            xs.clone().fold(f64::NEG_INFINITY, f64::max) - xs.fold(f64::INFINITY, f64::min)
        };
        ext(|q| q.0) + ext(|q| q.1)
    };
    (0..3)
        .map(around)
        .min_by(|a, b| span(a).total_cmp(&span(b)))
        .unwrap_or(w)
}

/// A triangle's parameter domain and its correspondence.
pub(crate) enum Piece {
    /// A parameter triangle mapped affinely onto the emitted one.
    Triangle([(f64, f64); 3], [Point3; 3]),
    /// A parameter rectangle mapped bilinearly onto four corners
    /// (`(u0, v0), (u1, v0), (u0, v1), (u1, v1)`).
    Rectangle((f64, f64), (f64, f64), [Point3; 4]),
}

impl Piece {
    /// The largest deviation on a 7 × 7 grid of the correspondence.
    pub(crate) fn sampled(&self, surface: &SurfaceGeometry) -> f64 {
        let n = 6;
        let mut worst = 0.0_f64;
        let mix = |c: &[Point3], w: &[f64]| {
            let mut out = Point3::new(0.0, 0.0, 0.0);
            for (p, k) in c.iter().zip(w) {
                out = Point3::new(out.x + k * p.x, out.y + k * p.y, out.z + k * p.z);
            }
            out
        };
        for i in 0..=n {
            for j in 0..=n {
                let (r, t) = (i as f64 / n as f64, j as f64 / n as f64);
                let (l, uv) = match self {
                    Piece::Triangle(w, p) => {
                        if r + t > 1.0 {
                            continue;
                        }
                        let k = [1.0 - r - t, r, t];
                        let uv = (
                            k[0] * w[0].0 + k[1] * w[1].0 + k[2] * w[2].0,
                            k[0] * w[0].1 + k[1] * w[1].1 + k[2] * w[2].1,
                        );
                        (mix(p, &k), uv)
                    }
                    Piece::Rectangle(u, v, c) => {
                        let k = [(1.0 - r) * (1.0 - t), r * (1.0 - t), (1.0 - r) * t, r * t];
                        (mix(c, &k), (u.0 + r * (u.1 - u.0), v.0 + t * (v.1 - v.0)))
                    }
                };
                worst = worst.max(l.distance_to(&surface.point_at(uv.0, uv.1)));
            }
        }
        worst
    }
}

/// The parameter pieces of a triangle. A vertex at a chart pole (a
/// sphere's pole, a cone's apex: the surface is the same point all
/// along one parameter line there) carries an arbitrary parameter, so a
/// pole fan's triangle `P A B` really stands for the rectangle from the
/// pole line to `A`'s and `B`'s line, mapped bilinearly with both pole
/// corners at `P`; neighbouring fan triangles then tile the pole's
/// neighbourhood. When `A` and `B` are off one line the domain is the
/// quad split in two triangles instead.
pub(crate) fn pole_pieces(
    surface: &SurfaceGeometry,
    w: [(f64, f64); 3],
    p: [Point3; 3],
    tolerance: f64,
) -> Vec<Piece> {
    // "The same point": a hundredth of the tolerance. The bound still
    // measures each piece's residual exactly (the corner
    // discrepancies); `near` only says the pole line is a point for
    // tiling, so a sliver the pieces leave is within it of the pole.
    let near = 1e-2 * tolerance;
    let at = |q: (f64, f64)| surface.point_at(q.0, q.1);
    for k in 0..3 {
        let (a, b) = ((k + 1) % 3, (k + 2) % 3);
        let pole = p[k];
        for axis in 0..2 {
            // The pole line runs along `axis` through w[k].
            let on_line = |other: (f64, f64)| {
                if axis == 0 {
                    (other.0, w[k].1)
                } else {
                    (w[k].0, other.1)
                }
            };
            let (pa, pb) = (on_line(w[a]), on_line(w[b]));
            if (pa == w[k] && pb == w[k])
                || at(pa).distance_to(&pole) > near
                || at(pb).distance_to(&pole) > near
            {
                continue;
            }
            let across = |q: (f64, f64)| if axis == 0 { q.1 } else { q.0 };
            if across(w[a]) == across(w[b]) {
                let (along_a, along_b) = if axis == 0 {
                    (w[a].0, w[b].0)
                } else {
                    (w[a].1, w[b].1)
                };
                let span = (along_a, along_b);
                let height = (across(w[k]), across(w[a]));
                return vec![if axis == 0 {
                    Piece::Rectangle(span, height, [pole, pole, p[a], p[b]])
                } else {
                    Piece::Rectangle(height, span, [pole, p[a], pole, p[b]])
                }];
            }
            return vec![
                Piece::Triangle([pa, w[a], w[b]], [pole, p[a], p[b]]),
                Piece::Triangle([pa, w[b], pb], [pole, p[b], pole]),
            ];
        }
    }
    vec![Piece::Triangle(w, p)]
}

/// C⁰ parameter lines and the periodic spline domain of `surface`.
pub(crate) fn c0_breaks(surface: &SurfaceGeometry) -> C0Breaks {
    let SurfaceGeometry::Nurbs(n) = surface else {
        // Analytic kinds are smooth; Gregory and polyhedral-spline
        // patches are single Bézier pieces in their exact forms.
        return C0Breaks::default();
    };
    let creases = |degree: usize, knots: &[f64], count: usize| -> Vec<f64> {
        let (lo, hi) = (knots[degree], knots[count]);
        let mut out = Vec::new();
        let mut k = degree + 1;
        while k < count {
            let x = knots[k];
            let mult = knots[k..count].iter().take_while(|&&y| y == x).count();
            if x > lo && x < hi && mult >= degree {
                out.push(x);
            }
            k += mult.max(1);
        }
        out
    };
    let mut breaks = C0Breaks {
        u: creases(n.degree_u, &n.knots_u, n.num_u),
        v: creases(n.degree_v, &n.knots_v, n.num_v),
        ..C0Breaks::default()
    };
    // A periodic chart's lifted triangles are cut at the domain's ends
    // and shifted into it.
    let domain = |degree: usize, knots: &[f64], p: Option<f64>| p.map(|p| (knots[degree], p));
    breaks.wrap_u = domain(n.degree_u, &n.knots_u, n.u_period());
    breaks.wrap_v = domain(n.degree_v, &n.knots_v, n.v_period());
    for (wrap, lines) in [
        (breaks.wrap_u, &mut breaks.u),
        (breaks.wrap_v, &mut breaks.v),
    ] {
        if let Some((start, p)) = wrap {
            lines.extend((-2..=3).map(|k| start + k as f64 * p));
        }
    }
    breaks
}
