//! Junction stations off the split vertex (PK-736,
//! `docs/architecture/fillet-open-tangent-chains.md`).
//!
//! A tangent-chain junction joins two strips on one exact section, and
//! the chain surgery keeps the junction's third edge (a crease, or a
//! periodic support's seam) by shortening it to its side's contact. At
//! the split vertex V that contact lies on the third edge only when the
//! section plane contains it. On an asymmetric tee (branch axis offset
//! from the main axis) the main cylinder's seam crosses the saddle
//! obliquely. Measured: the contact sits 0.317 off the seam at r = 2,
//! 0.870 at r = 6, and the surgery declined.
//!
//! The junction need not sit at V: the strips of one tangent chain meet
//! on whichever section is chosen. This module slides the section along
//! the chain (back along the arriving edge or on along the leaving one)
//! to the station where the contact lies ON the third edge. It solves by
//! bisection on the contact's signed lateral distance from the third
//! edge across its support. V then moves along the third edge to that
//! contact exactly as at a symmetric junction, so the topology is
//! unchanged.
//!
//! The search stays within 45% of each adjacent edge. No sign change
//! there keeps the original decline.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use super::tangent_chain_geometry::exact_section;

/// A section at a chain station: ball centre, contacts on P and Q, and
/// the section normal (the chain tangent there).
pub(super) struct Station {
    pub(super) center: Point3,
    pub(super) cp: Point3,
    pub(super) cq: Point3,
    pub(super) normal: Vec3,
}

/// The point and chain-direction tangent of `edge` at fraction `f` of
/// its span away from its end at `v`; `along` orients the tangent.
fn edge_station(
    brep: &BRep,
    edge: EdgeHandle,
    v: VertexHandle,
    f: f64,
    along: Vec3,
) -> Option<(Point3, Vec3)> {
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let vp = brep.vertices.get(v)?.position;
    let (t0, t1) = (curve.t_start, curve.t_end);
    let at_end = curve.geometry.point_at(t1).distance_to(&vp)
        <= curve.geometry.point_at(t0).distance_to(&vp);
    let (tv, tf) = if at_end { (t1, t0) } else { (t0, t1) };
    let t = tv + (tf - tv) * f;
    let tangent = curve.geometry.tangent_at(t).normalized()?;
    let tangent = if tangent.dot(&along) < 0.0 {
        tangent * -1.0
    } else {
        tangent
    };
    Some((curve.geometry.point_at(t), tangent))
}

/// The signed lateral distance of `c` from the curve of `third` across
/// `support` (positive on one fixed side), with the foot.
fn lateral(
    brep: &BRep,
    third: EdgeHandle,
    support: &SurfaceGeometry,
    c: Point3,
) -> Option<(f64, Point3)> {
    let curve = brep.edges.get(third)?.curve.as_ref()?;
    let s = curve.geometry.closest_point(&c)?;
    let foot = curve.geometry.point_at(s);
    let t = curve.geometry.tangent_at(s);
    let (u, v) = support.closest_point(&foot)?;
    let across = support.normal_at(u, v).cross(&t).normalized()?;
    Some(((c - foot).dot(&across), foot))
}

/// The station near V (arriving edge `ea`, leaving edge `eb`) whose
/// section contact on the third edge's side (`on_p`) lies on `third`,
/// within `band`.
#[allow(clippy::too_many_arguments)]
pub(super) fn station_on_third(
    brep: &BRep,
    (ea, eb): (EdgeHandle, EdgeHandle),
    v: VertexHandle,
    (sp, sq): (&SurfaceGeometry, &SurfaceGeometry),
    third: EdgeHandle,
    on_p: bool,
    radius: f64,
    (seed, along): (Point3, Vec3),
    band: f64,
) -> Option<Station> {
    let support = if on_p { sp } else { sq };
    // s < 0: back along the arriving edge; s > 0: on along the leaving.
    let at = |s: f64, seed: Point3| -> Option<(f64, Station)> {
        let (p, n) = if s < 0.0 {
            edge_station(brep, ea, v, -s, along)?
        } else {
            edge_station(brep, eb, v, s, along)?
        };
        let (center, cp, cq) = exact_section(sp, sq, p, n, radius, seed)?;
        let c = if on_p { cp } else { cq };
        let (d, _) = lateral(brep, third, support, c)?;
        Some((
            d,
            Station {
                center,
                cp,
                cq,
                normal: n,
            },
        ))
    };
    let (d0, s0) = at(0.0, seed)?;
    // Walk out both ways, nearest first, to the first sign change.
    const STEP: f64 = 0.01;
    let mut prev = [(0.0, d0, s0.center), (0.0, d0, s0.center)];
    let mut k = 1;
    while (k as f64) * STEP <= 0.45 + 1e-12 {
        let s = k as f64 * STEP;
        for (i, dir) in [1.0, -1.0].into_iter().enumerate() {
            let (lo, d_lo, c_lo) = prev[i];
            let Some((d, st)) = at(dir * s, c_lo) else {
                continue;
            };
            if d * d_lo <= 0.0 {
                let best = bisect(&at, dir, (lo, d_lo, c_lo), s, st)?;
                let c = if on_p { best.cp } else { best.cq };
                let (off, _) = lateral(brep, third, support, c)?;
                return (off.abs() <= band).then_some(best);
            }
            prev[i] = (s, d, st.center);
        }
        k += 1;
    }
    None
}

/// Bisect the station between `lo` (distance `d_lo`, ball `c_lo`) and
/// `hi` (whose section is `at_hi`) in direction `dir`.
fn bisect(
    at: &dyn Fn(f64, Point3) -> Option<(f64, Station)>,
    dir: f64,
    (lo, d_lo, c_lo): (f64, f64, Point3),
    hi: f64,
    at_hi: Station,
) -> Option<Station> {
    let (mut a, mut b, mut da) = (lo, hi, d_lo);
    let mut best = at_hi;
    let mut center = c_lo;
    for _ in 0..60 {
        let m = 0.5 * (a + b);
        let (dm, sm) = at(dir * m, center)?;
        center = sm.center;
        if dm * da <= 0.0 {
            b = m;
        } else {
            a = m;
            da = dm;
        }
        best = sm;
    }
    Some(best)
}
