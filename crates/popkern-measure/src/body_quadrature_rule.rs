//! The quadrature of one parameter triangle (`oriented_term`): Dunavant's
//! degree-5 rule, over the triangle clipped along a NURBS surface's knot
//! lines and refined adaptively where the area element varies strongly.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::chart::ChartDescriptor;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::Integrand;

/// A convex polygon in parameters.
type Polygon = Vec<(f64, f64)>;

/// Dunavant degree-5 rule: barycentric points and weights (sum 1).
const RULE: [([f64; 3], f64); 7] = [
    ([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], 0.225),
    (
        [
            0.059_715_871_789_77,
            0.470_142_064_105_115,
            0.470_142_064_105_115,
        ],
        0.132_394_152_788_506,
    ),
    (
        [
            0.470_142_064_105_115,
            0.059_715_871_789_77,
            0.470_142_064_105_115,
        ],
        0.132_394_152_788_506,
    ),
    (
        [
            0.470_142_064_105_115,
            0.470_142_064_105_115,
            0.059_715_871_789_77,
        ],
        0.132_394_152_788_506,
    ),
    (
        [
            0.797_426_985_353_087,
            0.101_286_507_323_456,
            0.101_286_507_323_456,
        ],
        0.125_939_180_544_827,
    ),
    (
        [
            0.101_286_507_323_456,
            0.797_426_985_353_087,
            0.101_286_507_323_456,
        ],
        0.125_939_180_544_827,
    ),
    (
        [
            0.101_286_507_323_456,
            0.101_286_507_323_456,
            0.797_426_985_353_087,
        ],
        0.125_939_180_544_827,
    ),
];

/// `kind`'s integrand over the oriented (u, v) triangle `uv`: the rule
/// over its pieces within one knot span (`knot_pieces`), each refined
/// adaptively where its area element varies strongly (`smooth_term`). The
/// atlas meshes a flat cap coarsely, its triangles spanning a quarter turn
/// of a revolved disc's rational chart: the rule alone read it 1e-3 short.
pub(in super::super) fn oriented_term(
    surface: &SurfaceGeometry,
    uv: [(f64, f64); 3],
    kind: Integrand,
) -> f64 {
    let pieces = knot_pieces(surface, uv);
    if pieces.is_empty() {
        return smooth_term(surface, uv, kind);
    }
    // Each convex piece fanned from its first corner, keeping the
    // triangle's orientation.
    pieces
        .iter()
        .flat_map(|poly| (1..poly.len() - 1).map(move |k| [poly[0], poly[k], poly[k + 1]]))
        .map(|t| smooth_term(surface, t, kind))
        .sum()
}

/// `uv` clipped along the interior knot lines of a NURBS surface it
/// crosses (both directions, the knots repeated by each evaluation
/// period), as convex pieces; empty when it crosses none. A NURBS
/// integrand is smooth only within a knot span (a rational circle's
/// quarter arcs meet with a kink in their speed): the rule over a
/// triangle spanning a knot read a revolved disc 3e-4 short, however
/// finely refined.
fn knot_pieces(surface: &SurfaceGeometry, uv: [(f64, f64); 3]) -> Vec<Polygon> {
    let SurfaceGeometry::Nurbs(nurbs) = surface else {
        return Vec::new();
    };
    let (pu, pv) = ChartDescriptor::evaluation_periods_of(surface);
    let lines = |knots: &[f64], period: Option<f64>, lo: f64, hi: f64| -> Vec<f64> {
        let (Some(&first), Some(&last)) = (knots.first(), knots.last()) else {
            return Vec::new();
        };
        let shifts: Vec<f64> = match period {
            Some(p) if p > 0.0 => {
                let k0 = ((lo - last) / p).floor() as i64;
                let k1 = ((hi - first) / p).ceil() as i64;
                (k0..=k1).map(|k| k as f64 * p).collect()
            }
            _ => vec![0.0],
        };
        let eps = 1e-12 * (hi - lo).abs().max(1e-300);
        let mut out: Vec<f64> = shifts
            .iter()
            .flat_map(|s| knots.iter().map(move |k| k + s))
            .filter(|&k| k > lo + eps && k < hi - eps)
            .collect();
        out.sort_by(f64::total_cmp);
        out.dedup();
        out
    };
    let (u_lo, u_hi) = (
        uv[0].0.min(uv[1].0).min(uv[2].0),
        uv[0].0.max(uv[1].0).max(uv[2].0),
    );
    let (v_lo, v_hi) = (
        uv[0].1.min(uv[1].1).min(uv[2].1),
        uv[0].1.max(uv[1].1).max(uv[2].1),
    );
    let cuts_u = lines(&nurbs.knots_u, pu, u_lo, u_hi);
    let cuts_v = lines(&nurbs.knots_v, pv, v_lo, v_hi);
    if cuts_u.is_empty() && cuts_v.is_empty() {
        return Vec::new();
    }
    let mut pieces = vec![uv.to_vec()];
    for (axis, cuts) in [(0usize, &cuts_u), (1usize, &cuts_v)] {
        for &c in cuts {
            pieces = pieces
                .into_iter()
                .flat_map(|poly| {
                    let (below, above) = clip(&poly, axis, c);
                    [below, above]
                })
                .filter(|poly| poly.len() >= 3)
                .collect();
        }
    }
    pieces
}

/// The convex polygon `poly` split by the line `coordinate[axis] = c`:
/// the part below and the part above, each in `poly`'s order.
fn clip(poly: &[(f64, f64)], axis: usize, c: f64) -> (Polygon, Polygon) {
    let at = |q: (f64, f64)| if axis == 0 { q.0 } else { q.1 };
    let (mut below, mut above) = (Vec::new(), Vec::new());
    let n = poly.len();
    for i in 0..n {
        let (p, q) = (poly[i], poly[(i + 1) % n]);
        let (dp, dq) = (at(p) - c, at(q) - c);
        if dp <= 0.0 {
            below.push(p);
        }
        if dp >= 0.0 {
            above.push(p);
        }
        if (dp < 0.0 && dq > 0.0) || (dp > 0.0 && dq < 0.0) {
            let t = dp / (dp - dq);
            let mut x = (p.0 + (q.0 - p.0) * t, p.1 + (q.1 - p.1) * t);
            if axis == 0 {
                x.0 = c;
            } else {
                x.1 = c;
            }
            below.push(x);
            above.push(x);
        }
    }
    (below, above)
}

/// `kind`'s integrand over a triangle within one knot span: the rule,
/// refined adaptively where the area element varies strongly.
fn smooth_term(surface: &SurfaceGeometry, uv: [(f64, f64); 3], kind: Integrand) -> f64 {
    let (whole, magnitude, (lo, hi)) = rule_term(surface, uv, kind);
    // Refined only where the area element varies across the triangle by
    // more than `ADAPTIVE_SPREAD` (read at the rule's own points, as near
    // a pole): elsewhere the rule's error is far under the agreement, and
    // checking every triangle's children made the measures five times
    // slower.
    if !hi.is_finite() || hi <= ADAPTIVE_SPREAD * lo {
        return whole;
    }
    refined_term(
        surface,
        uv,
        kind,
        (whole, ADAPTIVE_AGREEMENT * magnitude),
        ADAPTIVE_DEPTH,
    )
}

/// Ratio of the rule points' largest area element to their smallest past
/// which a triangle's rule value is refined.
const ADAPTIVE_SPREAD: f64 = 1.25;

/// Subdivision depth of the adaptive rule.
const ADAPTIVE_DEPTH: u32 = 4;

/// Agreement, relative to the integral of the integrand's magnitude over
/// the top triangle (an integral cancelling to nothing, a volume term
/// about a point on its surface, has no relative agreement), at which a
/// triangle's rule value stands.
const ADAPTIVE_AGREEMENT: f64 = 1e-12;

/// `whole` (the rule's value on `uv`) refined against its children
/// until they agree within `allowance`.
fn refined_term(
    surface: &SurfaceGeometry,
    uv: [(f64, f64); 3],
    kind: Integrand,
    (whole, allowance): (f64, f64),
    depth: u32,
) -> f64 {
    let mid = |a: (f64, f64), b: (f64, f64)| (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
    let (m01, m12, m20) = (mid(uv[0], uv[1]), mid(uv[1], uv[2]), mid(uv[2], uv[0]));
    let children = [
        [uv[0], m01, m20],
        [m01, uv[1], m12],
        [m20, m12, uv[2]],
        [m01, m12, m20],
    ];
    let parts = children.map(|c| rule_term(surface, c, kind).0);
    let sum: f64 = parts.iter().sum();
    if depth == 0 || (sum - whole).abs() <= allowance {
        return sum;
    }
    children
        .iter()
        .zip(parts)
        .map(|(c, part)| refined_term(surface, *c, kind, (part, 0.25 * allowance), depth - 1))
        .sum()
}

/// The degree-5 rule's value of `kind`'s integrand on `uv`, of its
/// magnitude, and the least and greatest area element at its points.
fn rule_term(
    surface: &SurfaceGeometry,
    uv: [(f64, f64); 3],
    kind: Integrand,
) -> (f64, f64, (f64, f64)) {
    let signed_area = 0.5
        * ((uv[1].0 - uv[0].0) * (uv[2].1 - uv[0].1) - (uv[2].0 - uv[0].0) * (uv[1].1 - uv[0].1));
    let (mut sum, mut magnitude) = (0.0, 0.0);
    let (mut lo, mut hi) = (f64::INFINITY, 0.0_f64);
    for (bary, w) in RULE {
        let u = bary[0] * uv[0].0 + bary[1] * uv[1].0 + bary[2] * uv[2].0;
        let v = bary[0] * uv[0].1 + bary[1] * uv[1].1 + bary[2] * uv[2].1;
        let n = surface.partial_u(u, v).cross(&surface.partial_v(u, v));
        let j = n.length();
        (lo, hi) = (lo.min(j), hi.max(j));
        let f = match kind {
            Integrand::Volume => {
                {
                    let p = surface.point_at(u, v);
                    Vec3::new(p.x, p.y, p.z)
                }
                .dot(&n)
                    / 3.0
            }
            Integrand::Area => n.length(),
        };
        sum += w * f;
        magnitude += w * f.abs();
    }
    (sum * signed_area, magnitude * signed_area.abs(), (lo, hi))
}
