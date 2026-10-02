//! Setback junctions between angularly adjacent strips (PK-527).
//!
//! Two adjacent strips share a support face F, and each has a rail on
//! it, its contact curve.
//! - **Convex face corner** (F's apex angle at the vertex below 180°):
//!   the rails cross inside both strips' ranges, and the setback
//!   corner sits at the crossing.
//! - **Reflex face corner** (a saddle): the rails diverge. Their
//!   virtual crossing lies behind the vertex, where a setback arc runs
//!   on a strip buried in material. The boundary instead gets a
//!   **bridge** on F, a curve from strip a's rail to strip b's rail.
//!
//! The classifier reads each rail's frame at its vertex end: the
//! direction `d` leaving the vertex, and the offset `m` from the edge
//! into F. At a convex corner each rail's direction points to the
//! other's interior side, so d_b·m_a > 0 and d_a·m_b > 0. Reflex and
//! near-flat corners (within [`FLAT_CORNER_SIN`]) take a bridge. Mixed
//! signs are not a face corner the model covers, and decline.
//!
//! The bridge leaves each strip at its rail point level with the
//! strip's other corner: the foot on the rail, whose normal plane holds
//! that corner. The strip's arc is then a pure cross-section.
//!
//! The bridge leaves p_a along rail a and arrives at p_b along rail b.
//! It is a Bézier whose control points run along each rail's line to
//! v, the crossing of the rails' tangent lines, so it follows each
//! rail through [`RAIL_FLAT_ORDER`] derivatives before turning near v.
//!
//! A cap corner needs compatible data from both sides. Each piece's
//! rim curvature must equal the other piece's cross-curvature.
//! - The strip's cross direction runs along its straight ruling, so
//!   the bridge must leave straight. An ellipse quarter (curvature 5.5
//!   at the notch) left a 0.118 crease. A quartic, straight through
//!   second order, left 3.9e-2. The rail-flat bridge leaves 8e-8.
//! - The strip's rim, its section, curves at 1/r, and the plane's own
//!   cross-curvature is 0. So the bridge's cross-curvature takes each
//!   strip's section curvature at its feet (`sides::bridge_cross_curvatures`).
//!   The cap stays G1 to the plane along the bridge, but not G2: G2 is
//!   unattainable at a reflex saddle's bridge feet.

use popkern_core::geometry::curve::{CurveEval, closest_point_on_curve};
use popkern_core::geometry::{NurbsCurve, Point3, SurfaceGeometry, Vec3};

use super::kp20_paper_border_corners::pair_is_coplanar;
use crate::fillet_surface::FilletSurfaceResult;

/// Order through which a bridge runs straight along each rail at its
/// feet. Its degree-2k Bézier has k + 1 control points on each rail's
/// line, so derivatives 1..k at each foot lie along the rail. That is
/// the corner compatibility a strip's straight ruling needs: a quartic
/// (k = 2) left the corner sectors creased at 3.9e-2. The KP20 rim
/// rows are degree 6, and k = 6 matches them.
const RAIL_FLAT_ORDER: usize = 6;

/// Fraction of each rail leg (foot to v) that the bridge's rail control
/// points span. The last leg into v is the remaining (1 − span): a
/// longer leg spreads the turn at the bridge's middle, the fan's spoke
/// point m.
///
/// Measured on the notched-cube saddles at k = 6, as the worst live
/// G1 sine, cap-cap / cap-boundary (PK-498; the `1 − cos θ` boundary
/// metric had hidden the boundary term):
///
/// | span | constant | variable |
/// |---|---|---|
/// | 0.7 | 1.2e-2 / 8.1e-3 | 1.4e-2 / 9.9e-3 |
/// | 0.8 | 8.6e-3 / 6.7e-3 | 1.0e-2 / 8.2e-3 |
/// | 0.85 | 7.1e-3 / 6.0e-3 | 8.5e-3 / 7.4e-3 |
/// | 0.9 | 5.7e-3 / 5.3e-3 | 6.9e-3 / 6.5e-3 |
/// | 0.95 | 7.6e-8 / 1.3e-1 | 5.5e-3 / 5.6e-3 |
/// | 1.0 | 8.2e-8 / 1.2e-1 | 4.2e-3 / 4.8e-3 |
///
/// Near a uniform span the constant saddle flips into a regime whose
/// cap-cap seams are exact but whose cap boundary creases 7° at a rim
/// piece corner by the spoke point. Assembly gates, corner turns and
/// section fits are unchanged across the flip, so its mechanism is
/// unexplained. 0.85 keeps a 0.1 margin from it. The fan commit
/// postcondition also gates the cap boundary, so a flip declines to
/// the ring rung instead of serving the crease.
const RAIL_FLAT_SPAN: f64 = 0.85;

/// `sin` of the apex-angle margin from 180° (5°) below which a face
/// corner counts as flat and takes a bridge.
const FLAT_CORNER_SIN: f64 = 0.087;

/// The kind of junction between two adjacent strips.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum JunctionKind {
    /// Convex face corner: the rails cross, and the setback corner sits
    /// at the crossing.
    Corner,
    /// Reflex or flat face corner: a bridge side on the shared face.
    Bridge,
}

/// A rail's frame at its vertex end.
struct RailFrame {
    /// Direction along the rail, leaving the vertex.
    away: Vec3,
    /// Offset from the edge into the shared face.
    inward: Vec3,
}

fn rail_frame(rail: &NurbsCurve, vertex: Point3) -> Option<RailFrame> {
    let (t0, t1) = rail.parameter_range();
    let (near, far) =
        if rail.point_at(t0).distance_to(&vertex) <= rail.point_at(t1).distance_to(&vertex) {
            (t0, t1)
        } else {
            (t1, t0)
        };
    let end = rail.point_at(near);
    let tangent = rail.derivative_at(near, 1).normalized()?;
    let away = if tangent.dot(&(rail.point_at(far) - end)) >= 0.0 {
        tangent
    } else {
        tangent * -1.0
    };
    let offset = end - vertex;
    let inward = (offset - away * offset.dot(&away)).normalized()?;
    Some(RailFrame { away, inward })
}

/// `|cos|` between a strip's normal at its rail and the rails' plane
/// normal above which the strip is tangent to that plane.
const TANGENT_PLANE_COS: f64 = 0.99;

/// The rails of `a` and `b` on their shared support face, as
/// `(rail_a, rail_b)`: the coplanar trim pair that both strips are
/// tangent to. Coplanarity alone is ambiguous. On a trihedral corner
/// the two rails on the NON-shared faces are coplanar too (both in
/// the plane at the contact offset, crossing on the edge), and taking
/// them read a convex cube corner as flat.
pub(super) fn shared_rails<'a>(
    a: &'a FilletSurfaceResult,
    b: &'a FilletSurfaceResult,
) -> Option<(&'a NurbsCurve, &'a NurbsCurve)> {
    [&a.trim_curve_1, &a.trim_curve_2]
        .into_iter()
        .flat_map(|ra| [(ra, &b.trim_curve_1), (ra, &b.trim_curve_2)])
        .find(|(ra, rb)| {
            pair_is_coplanar(ra, rb)
                && rails_plane_normal(ra, rb)
                    .is_some_and(|n| tangent_to(&a.surface, ra, n) && tangent_to(&b.surface, rb, n))
        })
}

fn midpoint_of(rail: &NurbsCurve) -> Point3 {
    let (t0, t1) = rail.parameter_range();
    rail.point_at(0.5 * (t0 + t1))
}

/// The unit normal of the plane two coplanar rails span.
fn rails_plane_normal(ra: &NurbsCurve, rb: &NurbsCurve) -> Option<Vec3> {
    let (t0, t1) = ra.parameter_range();
    let (a0, a1) = (ra.point_at(t0), ra.point_at(t1));
    let along = a1 - a0;
    let across = midpoint_of(rb) - a0;
    along.cross(&across).normalized().or_else(|| {
        let (s0, s1) = rb.parameter_range();
        along
            .cross(&(rb.point_at(s1) - rb.point_at(s0)))
            .normalized()
    })
}

/// Whether `strip` is tangent to the plane with normal `n` along its
/// `rail`: its normal at the rail's midpoint is parallel to `n`.
fn tangent_to(strip: &SurfaceGeometry, rail: &NurbsCurve, n: Vec3) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    let p = midpoint_of(rail);
    popkern_core::geometry::uv_projection::project_point_to_uv(strip, &p)
        .and_then(|(u, v)| strip.normal_at(u, v).normalized())
        .is_some_and(|normal| normal.dot(&n).abs() > TANGENT_PLANE_COS)
}

/// Classify the junction between adjacent strips `a` and `b` at
/// `vertex`.
pub(super) fn classify_junction(
    a: &FilletSurfaceResult,
    b: &FilletSurfaceResult,
    vertex: Point3,
) -> Result<JunctionKind, String> {
    let (rail_a, rail_b) =
        shared_rails(a, b).ok_or("adjacent strips have no coplanar rail pair")?;
    let (fa, fb) = (
        rail_frame(rail_a, vertex).ok_or("degenerate rail frame")?,
        rail_frame(rail_b, vertex).ok_or("degenerate rail frame")?,
    );
    let (sa, sb) = (fb.away.dot(&fa.inward), fa.away.dot(&fb.inward));
    if sa > FLAT_CORNER_SIN && sb > FLAT_CORNER_SIN {
        Ok(JunctionKind::Corner)
    } else if sa < FLAT_CORNER_SIN && sb < FLAT_CORNER_SIN {
        Ok(JunctionKind::Bridge)
    } else {
        Err(format!(
            "rails disagree on the face corner (sin {sa:.3} and {sb:.3})"
        ))
    }
}

/// A bridge side on a support face.
#[derive(Debug, Clone)]
pub struct BridgeGeometry {
    /// Where the bridge leaves strip a's rail.
    pub start: Point3,
    /// Where it meets strip b's rail.
    pub end: Point3,
    /// The bridge, `start` → `end`: a degree-12 Bézier with seven
    /// control points on each rail's line, meeting at the rails'
    /// crossing v (see [`RAIL_FLAT_ORDER`]). It runs along each rail
    /// through six derivatives, matching the strips' straight rulings
    /// at the corner.
    pub curve: NurbsCurve,
    /// The support plane the bridge lies in, spanned by the rails.
    pub plane: SurfaceGeometry,
}

/// How far a bowed bridge's end tangents turn towards its rails'
/// crossing: they point at the point this fraction of the way from the
/// chord's midpoint to the crossing (PK-619). At 1 each foot leaves
/// along its rail (a 90° cap corner), but on a steep apex the bridge
/// then pushes a long tongue into the cap. Measured on the 160-row
/// census with fractions 0.25 and 0.75 both in the ladder: 0.25 won
/// every row a bowed border served, 0.75 none.
const BRIDGE_BOW: f64 = 0.25;

/// Where a bowed bridge's inner controls sit along its end tangents, as
/// a fraction of the way to [`BRIDGE_BOW`]'s point (2/3 would be the
/// quadratic). Smaller turns the bridge nearer its feet and flattens its
/// middle, where the spokes meet the bridge: the quadratic's curved
/// middle forced the cap to twist at the spoke foot (pyramid-4 h3
/// variable G1 15.7 / r_min). Measured (honest reading, G1):
/// 2/3 -> 1/3 takes pyramid-4 h3 variable 15.7 -> 2.15, wide 7.30 ->
/// 2.70, pyramid-6 h0.4 wide 3.03 -> 1.03, pyramid-3 h1.2 variable
/// 3.84 -> 1.26; 0.2 declines on the steep trihedral apexes, 0.45 reads
/// 3.3 and 2.6 on the square ones.
const BRIDGE_TURN: f64 = 1.0 / 3.0;

/// A bridge's shape between its feet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum BridgeShape {
    /// Rail-flat through [`RAIL_FLAT_ORDER`], turning near the rails'
    /// crossing. At a reflex face corner the crossing lies beyond the
    /// cap, and the bridge is convex from the cap's side.
    RailFlat,
    /// The straight chord between the feet (PK-611). At a convex face
    /// corner a rail-flat bridge turns about the crossing, which lies
    /// inside the cap: the face's rounded corner pokes a tongue into
    /// the cap, and the fit's sectors fold at its tip (curvature
    /// 85–250 on the pentagon h0.8 wide apex, strips 7–28). The chord
    /// chamfers the face corner instead, so each foot is a convex cap
    /// corner of 90° plus the chord's angle to the rail.
    Chord,
    /// A cubic whose end tangents turn [`BRIDGE_BOW`] of the way towards
    /// the rails' crossing (PK-619), at the rail-flat bridge's degree:
    /// flat ([`BRIDGE_TURN`]) or round (the quadratic).
    Bowed(crate::types::BridgeBow),
}

/// The bridge on the shared face of strips `a` and `b`. It leaves `a`'s
/// rail level with `a`'s other corner `station_a`, and meets `b`'s rail
/// level with `station_b`.
pub(super) fn bridge_between(
    a: &FilletSurfaceResult,
    b: &FilletSurfaceResult,
    station_a: Point3,
    station_b: Point3,
    vertex: Point3,
    shape: BridgeShape,
) -> Result<BridgeGeometry, String> {
    let (rail_a, rail_b) =
        shared_rails(a, b).ok_or("adjacent strips have no coplanar rail pair")?;
    let foot = |rail: &NurbsCurve, station: Point3| -> Result<(Point3, Vec3), String> {
        let t = closest_point_on_curve(rail, &station).ok_or("rail foot not found")?;
        let tangent = rail
            .derivative_at(t, 1)
            .normalized()
            .ok_or("degenerate rail tangent")?;
        Ok((rail.point_at(t), tangent))
    };
    let (start, ta) = foot(rail_a, station_a)?;
    let (end, tb) = foot(rail_b, station_b)?;
    // The shoulder: the tangent lines' closest approach, which must be
    // an intersection (the rails are coplanar).
    let w = start - end;
    let (aa, ab, bb) = (ta.dot(&ta), ta.dot(&tb), tb.dot(&tb));
    let denom = aa * bb - ab * ab;
    if denom.abs() < 1.0e-12 {
        return Err("bridge rails are parallel".into());
    }
    let s = (ab * tb.dot(&w) - bb * ta.dot(&w)) / denom;
    let t = (aa * tb.dot(&w) - ab * ta.dot(&w)) / denom;
    let (on_a, on_b) = (start + ta * s, end + tb * t);
    let scale = start.distance_to(&end);
    if on_a.distance_to(&on_b) > 1.0e-6 * scale.max(1.0e-9) {
        return Err("bridge rails do not meet".into());
    }
    let shoulder = on_a + (on_b - on_a) * 0.5;
    // The shoulder must lie towards the vertex from both feet: past
    // each rail's vertex end, where the rails diverge.
    let towards_vertex = |p: Point3| (shoulder - p).dot(&(vertex - p)) > 0.0;
    if !(towards_vertex(start) && towards_vertex(end)) {
        return Err("bridge shoulder is not on the vertex side".into());
    }
    let normal = (shoulder - start)
        .cross(&(end - shoulder))
        .normalized()
        .ok_or("degenerate bridge shoulder")?;
    match shape {
        BridgeShape::Chord => return chord_bridge(start, end, normal),
        BridgeShape::Bowed(bow) => {
            let turn = match bow {
                crate::types::BridgeBow::Flat => BRIDGE_TURN,
                // A cubic with its inner controls two thirds of the way to
                // the quadratic's control is that quadratic.
                crate::types::BridgeBow::Round => 2.0 / 3.0,
            };
            // The bow's end tangents point at the control a quarter of
            // the way to the crossing; the inner controls sit a third
            // of the way along them, so the bridge turns near its feet
            // and runs nearly straight through its middle.
            let mid = start + (end - start) * 0.5;
            let towards = mid + (shoulder - mid) * BRIDGE_BOW;
            let controls = [
                start,
                start + (towards - start) * turn,
                end + (towards - end) * turn,
                end,
            ];
            return bowed_bridge(controls, normal);
        }
        BridgeShape::RailFlat => {}
    }
    let u_axis = (shoulder - start)
        .normalized()
        .ok_or("degenerate bridge leg")?;
    // Control points: RAIL_FLAT_ORDER + 1 of them along each rail,
    // meeting at the shoulder.
    let k = RAIL_FLAT_ORDER;
    let along = |p: Point3, i: usize| {
        if i == k {
            shoulder
        } else {
            p + (shoulder - p) * (RAIL_FLAT_SPAN * i as f64 / k as f64)
        }
    };
    let control_points: Vec<Point3> = (0..=k)
        .map(|i| along(start, i))
        .chain((0..k).rev().map(|i| along(end, i)))
        .collect();
    Ok(BridgeGeometry {
        start,
        end,
        curve: NurbsCurve {
            degree: 2 * k,
            knots: vec![0.0; 2 * k + 1]
                .into_iter()
                .chain(vec![1.0; 2 * k + 1])
                .collect(),
            weights: vec![1.0; control_points.len()],
            control_points,
        },
        plane: SurfaceGeometry::Plane {
            origin: start,
            normal,
            u_axis,
            v_axis: normal.cross(&u_axis),
        },
    })
}

/// The straight bridge `start` → `end` in the plane with unit
/// `normal`, at the rail-flat bridge's degree with uniform speed.
fn chord_bridge(start: Point3, end: Point3, normal: Vec3) -> Result<BridgeGeometry, String> {
    let u_axis = (end - start)
        .normalized()
        .ok_or("degenerate bridge chord")?;
    let degree = 2 * RAIL_FLAT_ORDER;
    let control_points: Vec<Point3> = (0..=degree)
        .map(|i| start + (end - start) * (i as f64 / degree as f64))
        .collect();
    Ok(BridgeGeometry {
        start,
        end,
        curve: NurbsCurve {
            degree,
            knots: vec![0.0; degree + 1]
                .into_iter()
                .chain(vec![1.0; degree + 1])
                .collect(),
            weights: vec![1.0; control_points.len()],
            control_points,
        },
        plane: SurfaceGeometry::Plane {
            origin: start,
            normal,
            u_axis,
            v_axis: normal.cross(&u_axis),
        },
    })
}

/// The cubic bridge with Bézier controls `c` in the plane with unit
/// `normal`, degree-elevated to the rail-flat bridge's degree.
fn bowed_bridge(c: [Point3; 4], normal: Vec3) -> Result<BridgeGeometry, String> {
    let (start, end) = (c[0], c[3]);
    let u_axis = (end - start)
        .normalized()
        .ok_or("degenerate bridge chord")?;
    let degree = 2 * RAIL_FLAT_ORDER;
    let binomial =
        |n: usize, k: usize| (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64);
    // Elevation from degree 3: control i of degree d is the binomial
    // blend of the cubic's four controls.
    let control_points: Vec<Point3> = (0..=degree)
        .map(|i| {
            let mut p = Vec3::new(0.0, 0.0, 0.0);
            for (j, q) in c.into_iter().enumerate() {
                if j <= i && i - j <= degree - 3 {
                    let w = binomial(3, j) * binomial(degree - 3, i - j) / binomial(degree, i);
                    p = p + Vec3::new(q.x, q.y, q.z) * w;
                }
            }
            Point3::new(p.x, p.y, p.z)
        })
        .collect();
    Ok(BridgeGeometry {
        start,
        end,
        curve: NurbsCurve {
            degree,
            knots: vec![0.0; degree + 1]
                .into_iter()
                .chain(vec![1.0; degree + 1])
                .collect(),
            weights: vec![1.0; control_points.len()],
            control_points,
        },
        plane: SurfaceGeometry::Plane {
            origin: start,
            normal,
            u_axis,
            v_axis: normal.cross(&u_axis),
        },
    })
}

#[cfg(test)]
#[path = "kp20_paper_border_junctions_tests.rs"]
mod tests;
