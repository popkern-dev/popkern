//! A NURBS face as the Bézier patches its trim reaches, and two faces
//! proved apart by subdividing the patches still near (PK-1155).
//!
//! A NURBS face's grid box can fall short of the face, so no closed
//! form clears a pair with one, and the shell search, which projects
//! onto the surface from every seed, cost seconds a pair on the
//! cycloidal disc's walls: two helical flanks a crest 0.05 wide apart,
//! and pieces of one closed wall whose whole surfaces coincide, so only
//! their trims are apart. A rational patch of positive weights lies in
//! the hull of its control points, and a face on a surface whose
//! parameters form a rectangle lies inside the uv box of its boundary
//! loops' pcurves (`pcurve_box`: closed form), so the patches of the
//! knot spans meeting that box hold the face. Where the box runs past
//! the surface's domain (a pcurve lifted across a seam) that direction
//! is not restricted. Each span's patch comes from its own p + 1
//! control points and 2p + 2 knots, its end knots raised to
//! multiplicity p (Boehm's insertion, in homogeneous coordinates), and
//! is cut to the part of the span inside the box (de Casteljau).
//!
//! Two hulls are apart when some axis separates their points by more
//! than the band: the box axes, an axis across either patch (the cross
//! product of its net's diagonal sums), their sum and difference (two
//! flanks either side of a crest), and the line between the patches'
//! middles. A pair no axis separates halves the larger patch along its
//! longer direction (de Casteljau; the hulls close on the surface
//! quadratically) down to a depth, and a pair left undecided there is
//! for the search.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::parallel::ff_parallel::ff_parallel_reach::{box_clear_of_face, box_gap};
use crate::pcurve_box::{UvBox, contour_uv_box};

/// Halvings of a pair's patches before it is left to the search.
const DEPTH: usize = 10;

/// A rational Bézier patch: `(nu) x (nv)` homogeneous control points
/// `(w x, w y, w z, w)`, u-major, with the box, middle and an axis
/// across of its points.
#[derive(Clone, Debug)]
pub(super) struct Patch {
    nu: usize,
    nv: usize,
    net: Vec<[f64; 4]>,
    points: Vec<Point3>,
    /// The box of the control points: it holds the patch.
    pub(super) bx: Aabb,
    middle: Point3,
    across: Option<Vec3>,
}

impl Patch {
    fn new(nu: usize, nv: usize, net: Vec<[f64; 4]>) -> Option<Self> {
        if net.len() != nu * nv || nu < 2 || nv < 2 || net.iter().any(|h| h[3] <= 0.0) {
            return None;
        }
        let points: Vec<Point3> = net
            .iter()
            .map(|h| Point3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]))
            .collect();
        let at = |a: usize, b: usize| points[a * nv + b];
        let (p, q) = (nu - 1, nv - 1);
        let du = (at(p, 0) - at(0, 0)) + (at(p, q) - at(0, q));
        let dv = (at(0, q) - at(0, 0)) + (at(p, q) - at(p, 0));
        let base = points[0];
        let sum = points
            .iter()
            .fold(Vec3::new(0.0, 0.0, 0.0), |acc, x| acc + (*x - base));
        Some(Self {
            nu,
            nv,
            bx: Aabb::from_points(&points)?,
            middle: base + sum * (1.0 / points.len() as f64),
            across: du.cross(&dv).normalized(),
            points,
            net,
        })
    }

    /// The patch's control polygon's longest extent along u and along v.
    fn extents(&self) -> (f64, f64) {
        let at = |a: usize, b: usize| self.points[a * self.nv + b];
        let along_u = (0..self.nv)
            .map(|b| at(0, b).distance_to(&at(self.nu - 1, b)))
            .fold(0.0, f64::max);
        let along_v = (0..self.nu)
            .map(|a| at(a, 0).distance_to(&at(a, self.nv - 1)))
            .fold(0.0, f64::max);
        (along_u, along_v)
    }

    /// The two halves at the parameter middle, along u when `in_u`.
    fn halves(&self, in_u: bool) -> Option<(Patch, Patch)> {
        self.split(in_u, 0.5)
    }

    /// The part over `[t0, t1]` of its parameter range, along u when
    /// `in_u` (itself when that is the whole range).
    fn clipped(self, in_u: bool, (t0, t1): (f64, f64)) -> Option<Patch> {
        let below = if t1 < 1.0 {
            self.split(in_u, t1)?.0
        } else {
            self
        };
        if t0 > 0.0 && t1 > 0.0 {
            Some(below.split(in_u, t0 / t1)?.1)
        } else {
            Some(below)
        }
    }

    /// The two parts either side of the parameter `t` (de Casteljau),
    /// along u when `in_u`.
    fn split(&self, in_u: bool, t: f64) -> Option<(Patch, Patch)> {
        let (nu, nv) = (self.nu, self.nv);
        let (mut lo, mut hi) = (self.net.clone(), self.net.clone());
        let lines = if in_u { nv } else { nu };
        let len = if in_u { nu } else { nv };
        let index = |line: usize, k: usize| if in_u { k * nv + line } else { line * nv + k };
        for line in 0..lines {
            let mut row: Vec<[f64; 4]> = (0..len).map(|k| self.net[index(line, k)]).collect();
            for k in 0..len {
                lo[index(line, k)] = row[0];
                hi[index(line, len - 1 - k)] = row[len - 1 - k];
                for m in 0..len - 1 - k {
                    row[m] = std::array::from_fn(|c| row[m][c] + t * (row[m + 1][c] - row[m][c]));
                }
            }
        }
        Some((Patch::new(nu, nv, lo)?, Patch::new(nu, nv, hi)?))
    }
}

/// The Bézier patches holding `face` (on `s`): `None` when a weight is
/// not positive, a boundary fin has no pcurve with a closed-form box, or
/// the surface's net is malformed.
pub(super) fn face_patches(brep: &BRep, face: FaceHandle, s: &NurbsSurface) -> Option<Vec<Patch>> {
    window_patches(s, &trim_window(brep, face, s)?)
}

/// The Bézier patches of the knot spans of `s` meeting the parameter
/// window `((u0, u1), (v0, v1))`, u-major in span order.
fn window_patches(s: &NurbsSurface, window: &((f64, f64), (f64, f64))) -> Option<Vec<Patch>> {
    let (p, q) = (s.degree_u, s.degree_v);
    if p == 0
        || q == 0
        || s.num_u <= p
        || s.num_v <= q
        || s.knots_u.len() != s.num_u + p + 1
        || s.knots_v.len() != s.num_v + q + 1
        || s.control_points.len() != s.num_u * s.num_v
        || s.weights.len() != s.control_points.len()
        || s.weights.iter().any(|&w| w <= 0.0)
    {
        return None;
    }
    let meets = |knots: &[f64], i: usize, (lo, hi): (f64, f64)| {
        knots[i] < knots[i + 1] && knots[i + 1] >= lo && knots[i] <= hi
    };
    let mut out = Vec::new();
    for i in (p..s.num_u).filter(|&i| meets(&s.knots_u, i, window.0)) {
        for j in (q..s.num_v).filter(|&j| meets(&s.knots_v, j, window.1)) {
            let h = |a: usize, b: usize| {
                let k = a * s.num_v + b;
                let (c, w) = (s.control_points[k], s.weights[k]);
                [c.x * w, c.y * w, c.z * w, w]
            };
            // Each of the q + 1 columns to Bézier in u, then each row in v.
            let mut cols: Vec<Vec<[f64; 4]>> = Vec::with_capacity(q + 1);
            for b in j - q..=j {
                let col: Vec<[f64; 4]> = (i - p..=i).map(|a| h(a, b)).collect();
                cols.push(span_bezier(&col, &s.knots_u[i - p..=i + p + 1], p)?);
            }
            let mut net = vec![[0.0; 4]; (p + 1) * (q + 1)];
            for a in 0..=p {
                let row: Vec<[f64; 4]> = cols.iter().map(|c| c[a]).collect();
                let bez = span_bezier(&row, &s.knots_v[j - q..=j + q + 1], q)?;
                for (b, x) in bez.into_iter().enumerate() {
                    net[a * (q + 1) + b] = x;
                }
            }
            // Only the part of the span inside the window: a face cut
            // across a span holds that part alone (a thread flank either
            // side of a slot through one span, each face holding the
            // whole span, was never apart from the other, PK-1172).
            let part = |knots: &[f64], k: usize, (lo, hi): (f64, f64)| {
                let (a, b) = (knots[k], knots[k + 1]);
                (
                    ((lo - a) / (b - a)).clamp(0.0, 1.0),
                    ((hi - a) / (b - a)).clamp(0.0, 1.0),
                )
            };
            let (tu, tv) = (part(&s.knots_u, i, window.0), part(&s.knots_v, j, window.1));
            if tu.1 <= tu.0 || tv.1 <= tv.0 {
                continue;
            }
            let patch = Patch::new(p + 1, q + 1, net)?
                .clipped(true, tu)?
                .clipped(false, tv)?;
            out.push(patch);
        }
    }
    (!out.is_empty()).then_some(out)
}

/// The Bézier points of the middle span of a degree-`p` curve given by
/// its `p + 1` homogeneous points and `2p + 2` local knots: the span's
/// end knots (`knots[p]`, `knots[p + 1]`) raised to multiplicity `p`.
fn span_bezier(points: &[[f64; 4]], knots: &[f64], p: usize) -> Option<Vec<[f64; 4]>> {
    let (a, b) = (knots[p], knots[p + 1]);
    if points.len() != p + 1 || knots.len() != 2 * p + 2 || a >= b {
        return None;
    }
    let (mut pts, mut u) = (points.to_vec(), knots.to_vec());
    for t in [a, b] {
        while u.iter().filter(|&&x| x == t).count() < p {
            // The span of the local domain [a, b] at t: U[k] <= a < U[k + 1]
            // at its start, U[k] < b <= U[k + 1] at its end (the half-open
            // span past b holds no p + 1 points).
            let k = if t == a {
                (0..u.len() - 1)
                    .rev()
                    .find(|&k| u[k] <= t && t < u[k + 1])?
            } else {
                (0..u.len() - 1).find(|&k| u[k] < t && t <= u[k + 1])?
            };
            if k < p || k >= pts.len() {
                return None;
            }
            let mut next = Vec::with_capacity(pts.len() + 1);
            for jj in 0..=pts.len() {
                next.push(if jj + p <= k {
                    pts[jj]
                } else if jj > k {
                    pts[jj - 1]
                } else {
                    let den = u[jj + p] - u[jj];
                    let al = if den > 0.0 { (t - u[jj]) / den } else { 0.0 };
                    std::array::from_fn(|c| (1.0 - al) * pts[jj - 1][c] + al * pts[jj][c])
                });
            }
            pts = next;
            u.insert(k + 1, t);
        }
    }
    let s = (0..u.len() - 1).find(|&s| u[s] == a && u[s + 1] == b)?;
    (s >= p && s < pts.len()).then(|| pts[s - p..=s].to_vec())
}

/// The piece of `s` over the knot spans `face`'s trim reaches, with the
/// same parameters (`split_at_u` / `split_at_v` at span ends: exact),
/// or `None` when the face reaches every span or its window is unknown.
pub(super) fn face_piece(brep: &BRep, face: FaceHandle, s: &NurbsSurface) -> Option<NurbsSurface> {
    piece_over(s, trim_window(brep, face, s)?)
}

/// The piece of `s` over the knot spans the window meets, `None` when it
/// meets every span.
fn piece_over(
    s: &NurbsSurface,
    ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64)),
) -> Option<NurbsSurface> {
    // The window widened to the knots about it.
    let widen = |knots: &[f64], degree: usize, count: usize, lo: f64, hi: f64| {
        let (d0, d1) = (knots[degree], knots[count]);
        let below = knots[degree..=count]
            .iter()
            .copied()
            .filter(|&k| k <= lo)
            .fold(d0, f64::max);
        let above = knots[degree..=count]
            .iter()
            .copied()
            .filter(|&k| k >= hi)
            .fold(d1, f64::min);
        ((below > d0).then_some(below), (above < d1).then_some(above))
    };
    let (us, vs) = (
        widen(&s.knots_u, s.degree_u, s.num_u, u0, u1),
        widen(&s.knots_v, s.degree_v, s.num_v, v0, v1),
    );
    if us == (None, None) && vs == (None, None) {
        return None;
    }
    let mut piece = s.clone();
    if let Some(a) = us.0 {
        piece = piece.split_at_u(a).ok()?.1;
    }
    if let Some(b) = us.1 {
        piece = piece.split_at_u(b).ok()?.0;
    }
    if let Some(a) = vs.0 {
        piece = piece.split_at_v(a).ok()?.1;
    }
    if let Some(b) = vs.1 {
        piece = piece.split_at_v(b).ok()?.0;
    }
    Some(piece)
}

/// The face's parameter window `((u0, u1), (v0, v1))`: its boundary's
/// uv box, each direction the whole domain where that box runs past it.
fn trim_window(
    brep: &BRep,
    face: FaceHandle,
    s: &NurbsSurface,
) -> Option<((f64, f64), (f64, f64))> {
    let f = brep.faces.get(face)?;
    let mut pcurves = Vec::new();
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let lp = brep.loops.get(l)?;
        if lp.vertex.is_some() {
            continue;
        }
        let start = lp.half_edge;
        let mut h = start;
        for _ in 0..=brep.half_edges.len() {
            let fin = brep.half_edges.get(h)?;
            pcurves.push(fin.pcurve.as_ref()?);
            h = fin.next;
            if h == start {
                break;
            }
        }
    }
    let uv: UvBox = contour_uv_box(pcurves)?;
    let domain = |knots: &[f64], degree: usize, count: usize| (knots[degree], knots[count]);
    let (du, dv) = (
        domain(&s.knots_u, s.degree_u, s.num_u),
        domain(&s.knots_v, s.degree_v, s.num_v),
    );
    let within = |lo: f64, hi: f64, (d0, d1): (f64, f64)| {
        let slack = 1e-9 * (d1 - d0).abs().max(1.0);
        if lo >= d0 - slack && hi <= d1 + slack {
            (lo.max(d0), hi.min(d1))
        } else {
            (d0, d1)
        }
    };
    Some((within(uv.0.0, uv.1.0, du), within(uv.0.1, uv.1.1, dv)))
}

/// The box of the patches' boxes: a box holding the face.
pub(super) fn hull(patches: &[Patch]) -> Option<Aabb> {
    let corners: Vec<_> = patches.iter().flat_map(|s| [s.bx.min, s.bx.max]).collect();
    Aabb::from_points(&corners)
}

/// Whether the ball at a point of a radius lies off a face.
pub(super) type OffFace<'a> = dyn Fn(Point3, f64) -> bool + 'a;

/// The other face of a pair, as the patch test reads it.
pub(super) enum Other<'a> {
    /// A face on `surface` with its box (trusted only when the surface's
    /// boxes are known to hold their faces) and, for a planar face, a
    /// test of whether the ball about a point of its plane's band, of a
    /// radius, lies off the face itself.
    Face(&'a SurfaceGeometry, &'a Aabb, Option<&'a OffFace<'a>>),
    /// A NURBS face's patches.
    Patches(&'a [Patch]),
}

/// Whether a NURBS face with patches `mine` provably lies farther than
/// `pad` from `other`.
pub(super) fn patches_apart(mine: &[Patch], other: &Other, pad: f64) -> bool {
    mine.iter().all(|a| match other {
        Other::Face(surface, face_box, off) => {
            patch_clear_of_face(a, (surface, face_box, *off), pad, DEPTH)
        }
        Other::Patches(theirs) => theirs.iter().all(|b| pair_apart(a, b, pad, DEPTH)),
    })
}

fn patch_clear_of_face(
    a: &Patch,
    face: (&SurfaceGeometry, &Aabb, Option<&OffFace>),
    pad: f64,
    depth: usize,
) -> bool {
    let (surface, face_box, off) = face;
    if box_clear_of_face(&a.bx, surface, face_box, pad) {
        return true;
    }
    // The box's ball (about its middle, of its half-diagonal) off a
    // planar face by the pad.
    let middle = a.bx.min + (a.bx.max - a.bx.min) * 0.5;
    if off.is_some_and(|off| off(middle, pad + 0.5 * a.bx.min.distance_to(&a.bx.max))) {
        return true;
    }
    if depth == 0 {
        return false;
    }
    let (eu, ev) = a.extents();
    a.halves(eu >= ev).is_some_and(|(x, y)| {
        patch_clear_of_face(&x, face, pad, depth - 1)
            && patch_clear_of_face(&y, face, pad, depth - 1)
    })
}

fn pair_apart(a: &Patch, b: &Patch, pad: f64, depth: usize) -> bool {
    if hulls_apart(a, b, pad) {
        return true;
    }
    if depth == 0 {
        return false;
    }
    let size = |p: &Patch| p.bx.min.distance_to(&p.bx.max);
    let (big, small) = if size(a) >= size(b) { (a, b) } else { (b, a) };
    let (eu, ev) = big.extents();
    big.halves(eu >= ev).is_some_and(|(x, y)| {
        pair_apart(&x, small, pad, depth - 1) && pair_apart(&y, small, pad, depth - 1)
    })
}

/// Whether the control hulls of two patches lie farther than `pad`
/// apart along some axis (any axis separating the point sets by more
/// than `pad` separates their hulls by as much).
fn hulls_apart(a: &Patch, b: &Patch, pad: f64) -> bool {
    if box_gap(&a.bx, &b.bx) > pad {
        return true;
    }
    let sum = a
        .across
        .zip(b.across)
        .and_then(|(x, y)| (x + y).normalized());
    let diff = a
        .across
        .zip(b.across)
        .and_then(|(x, y)| (x - y).normalized());
    let axes = [
        a.across,
        b.across,
        sum,
        diff,
        (b.middle - a.middle).normalized(),
    ];
    axes.into_iter().flatten().any(|d| {
        let range = |s: &Patch| {
            s.points.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| {
                let x = Vec3::new(p.x, p.y, p.z).dot(&d);
                (lo.min(x), hi.max(x))
            })
        };
        let ((a0, a1), (b0, b1)) = (range(a), range(b));
        b0 - a1 > pad || a0 - b1 > pad
    })
}

#[path = "self_contact_patches_curve.rs"]
mod curve;
pub(super) use curve::curve_apart;

#[cfg(test)]
#[path = "self_contact_patches_tests.rs"]
mod tests;
