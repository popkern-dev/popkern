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
//!   **bridge** on F, a conic from strip a's rail to strip b's rail.
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
//! The bridge is a quartic that leaves p_a along rail a and arrives
//! at p_b along rail b, both with zero curvature. Its middle control
//! point is v, the crossing of the rails' tangent lines, and the
//! control points either side of v are halfway along each rail.
//!
//! A cap corner needs compatible data from both sides. Each piece's
//! rim curvature must equal the other piece's cross-curvature.
//! - The strip's cross direction runs along its straight ruling, so
//!   the bridge must leave with zero curvature: an ellipse quarter
//!   (curvature 5.5 at the notch) kept one corner crease.
//! - The strip's rim, its section, curves at 1/r, and the plane's own
//!   cross-curvature is 0. So the bridge's cross-curvature takes each
//!   strip's section curvature at its feet ([`bridge_cross_curvatures`]).
//!   The cap stays G1 to the plane along the bridge, but not G2: G2 is
//!   unattainable at a reflex saddle's bridge feet.

use popkern_core::geometry::curve::{CurveEval, closest_point_on_curve};
use popkern_core::geometry::{NurbsCurve, Point3, SurfaceGeometry, Vec3};

use super::kp20_paper_border_corners::pair_is_coplanar;
use crate::fillet_surface::FilletSurfaceResult;

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
    /// The bridge, `start` → `end`: the quartic with control points
    /// start, ½(start + v), v, ½(end + v), end. Collinear triples at
    /// both ends leave each rail along it with zero curvature,
    /// matching the strips' straight rulings at the corner.
    pub curve: NurbsCurve,
    /// The support plane the bridge lies in, spanned by the rails.
    pub plane: SurfaceGeometry,
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
    let u_axis = (shoulder - start)
        .normalized()
        .ok_or("degenerate bridge leg")?;
    let half = |p: Point3| p + (shoulder - p) * 0.5;
    Ok(BridgeGeometry {
        start,
        end,
        curve: NurbsCurve {
            degree: 4,
            knots: vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![start, half(start), shoulder, half(end), end],
            weights: vec![1.0; 5],
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

/// Which junctions a setback border may take.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SetbackJunctionPolicy {
    /// Every junction is a corner (the convex setback; reflex rails
    /// then meet at a gapped corner).
    CornersOnly,
    /// Reflex face corners take a bridge side (the saddle setback).
    BridgeReflexFaces,
}

/// One side of the setback boundary.
#[derive(Debug, Clone)]
pub(super) enum SideKind {
    /// A strip's section path.
    Strip {
        /// Index into the fillet surfaces.
        fillet_index: usize,
    },
    /// A bridge on a reflex support face.
    Bridge(Box<BridgeGeometry>),
}

/// One side of the setback boundary, between two junction points.
#[derive(Debug, Clone)]
pub(super) struct SetbackSide {
    pub kind: SideKind,
    pub start: Point3,
    pub end: Point3,
}

/// The setback boundary's sides in angular order. Strip `order[k]`
/// runs from junction k−1 to junction k. A bridge junction inserts a
/// bridge side between its two strips, and each strip then ends at the
/// bridge's foot on its rail. `corners[k]` is junction k's rail
/// crossing, used at corner junctions and as the anchoring station of
/// a bridge.
pub(super) fn setback_sides(
    fillet_surfaces: &[FilletSurfaceResult],
    order: &[usize],
    corners: &[(f64, Point3)],
    vertex: Point3,
    policy: SetbackJunctionPolicy,
) -> Result<Vec<SetbackSide>, String> {
    let n = order.len();
    let kinds: Vec<JunctionKind> = match policy {
        SetbackJunctionPolicy::CornersOnly => vec![JunctionKind::Corner; n],
        SetbackJunctionPolicy::BridgeReflexFaces => (0..n)
            .map(|k| {
                classify_junction(
                    &fillet_surfaces[order[k]],
                    &fillet_surfaces[order[(k + 1) % n]],
                    vertex,
                )
                // An unclassifiable junction keeps the corner.
                .unwrap_or(JunctionKind::Corner)
            })
            .collect(),
    };
    let bridges: Vec<Option<BridgeGeometry>> = (0..n)
        .map(|k| {
            if kinds[k] != JunctionKind::Bridge {
                return Ok(None);
            }
            let (prev, next) = ((k + n - 1) % n, (k + 1) % n);
            if kinds[prev] == JunctionKind::Bridge || kinds[next] == JunctionKind::Bridge {
                return Err("a strip between two reflex faces has no anchoring corner".to_string());
            }
            bridge_between(
                &fillet_surfaces[order[k]],
                &fillet_surfaces[order[next]],
                corners[prev].1,
                corners[next].1,
                vertex,
            )
            .map(Some)
        })
        .collect::<Result<_, _>>()?;
    let mut sides = Vec::with_capacity(n + bridges.iter().flatten().count());
    for k in 0..n {
        let prev = (k + n - 1) % n;
        let start = match &bridges[prev] {
            Some(bridge) => bridge.end,
            None => corners[prev].1,
        };
        let end = match &bridges[k] {
            Some(bridge) => bridge.start,
            None => corners[k].1,
        };
        sides.push(SetbackSide {
            kind: SideKind::Strip {
                fillet_index: order[k],
            },
            start,
            end,
        });
        if let Some(bridge) = &bridges[k] {
            sides.push(SetbackSide {
                kind: SideKind::Bridge(Box::new(bridge.clone())),
                start: bridge.start,
                end: bridge.end,
            });
        }
    }
    Ok(sides)
}

/// A side's path on its own surface.
#[derive(Debug, Clone)]
pub(super) enum SidePath {
    /// A strip's section path.
    Strip(super::section_path::SectionPath),
    /// A bridge quartic on its support plane.
    Bridge(Box<BridgeGeometry>),
}

impl SidePath {
    /// uv at path parameter `t ∈ [0, 1]` on [`Self::surface`].
    pub(super) fn uv(&self, t: f64) -> (f64, f64) {
        match self {
            Self::Strip(path) => path.uv(t),
            Self::Bridge(bridge) => plane_uv(&bridge.plane, bridge.curve.point_at(t)),
        }
    }

    /// The surface the side lies on.
    pub(super) fn surface<'a>(
        &'a self,
        fillet_surfaces: &'a [FilletSurfaceResult],
    ) -> &'a SurfaceGeometry {
        match self {
            Self::Strip(path) => &fillet_surfaces[path.fillet_index].surface,
            Self::Bridge(bridge) => &bridge.plane,
        }
    }

    /// The strip the side runs on; `None` for a bridge.
    pub(super) fn fillet_index(&self) -> Option<usize> {
        match self {
            Self::Strip(path) => Some(path.fillet_index),
            Self::Bridge(_) => None,
        }
    }
}

fn plane_uv(plane: &SurfaceGeometry, p: Point3) -> (f64, f64) {
    match plane {
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => ((p - *origin).dot(u_axis), (p - *origin).dot(v_axis)),
        _ => (0.0, 0.0),
    }
}

/// `uv` moved onto the periodic branch nearest `anchor` in each
/// periodic direction of `surface`. A strip spans under a half turn,
/// so its arc is the short way round. The raw corner projections
/// straddled the cylinder seam on the saddle's convex strips, and the
/// path ran 270° through material.
pub(super) fn unwrap_onto(
    surface: &SurfaceGeometry,
    anchor: (f64, f64),
    uv: (f64, f64),
) -> (f64, f64) {
    use popkern_core::geometry::surface::SurfaceEval;
    let (period_u, period_v) = surface.parameter_periodicity();
    let near = |value: f64, anchor: f64, period: Option<f64>| match period {
        Some(t) if t > 0.0 => value - ((value - anchor) / t).round() * t,
        _ => value,
    };
    (
        near(uv.0, anchor.0, period_u),
        near(uv.1, anchor.1, period_v),
    )
}

/// Give each bridge side's cross-curvature samples the neighbouring
/// strips' section curvature at its feet, blended along the bridge by
/// the same smootherstep the march uses. The plane's own
/// cross-curvature is 0, and without this the corner blocks are
/// incompatible. The normal part follows the march convention
/// `n·d2 = 1.2·κ·|jet|²`, and the in-plane part (row shaping) is kept.
pub(super) fn bridge_cross_curvatures(
    paths: &[SidePath],
    arcs: &[super::Kp20SetbackArc],
    shims: &[super::OrientedSegment],
    samples: &mut [[super::CrossFieldSamples; 2]],
) {
    let n = paths.len();
    for k in 0..n {
        if matches!(paths[k], SidePath::Bridge(_)) {
            let (prev, next) = ((k + n - 1) % n, (k + 1) % n);
            bridge_side_cross_curvatures(
                &shims[k],
                [&shims[prev], &shims[next]],
                [arcs[k].corner_prev, arcs[k].corner_next],
                &mut samples[k],
            );
        }
    }
}

/// [`bridge_cross_curvatures`] for one bridge side.
fn bridge_side_cross_curvatures(
    bridge: &super::OrientedSegment,
    strips: [&super::OrientedSegment; 2],
    feet: [Point3; 2],
    samples: &mut [super::CrossFieldSamples; 2],
) {
    let Some(normal) = bridge.surface_normal(feet[0]).and_then(|n| n.normalized()) else {
        return;
    };
    let end_jet = |half: usize| -> Option<Vec3> {
        let jets = &samples[half].jets;
        if half == 0 { jets.first() } else { jets.last() }.map(|(_, jet)| *jet)
    };
    let kappa = |end: usize| -> f64 {
        end_jet(end)
            .and_then(|jet| strips[end].normal_curvature(feet[end], jet * -1.0, normal))
            .unwrap_or(0.0)
    };
    let (kappa_start, kappa_end) = (kappa(0), kappa(1));
    for (half, piece) in samples.iter_mut().enumerate() {
        for (i, (t, sample)) in piece.curvatures.iter_mut().enumerate() {
            let Some(&(_, jet)) = piece.jets.get(i) else {
                continue;
            };
            let s = 0.5 * (half as f64 + *t);
            let w = s * s * s * (10.0 + s * (-15.0 + 6.0 * s));
            let k = kappa_start * (1.0 - w) + kappa_end * w;
            let in_plane = *sample - normal * sample.dot(&normal);
            *sample = in_plane + normal * (1.2 * k * jet.dot(&jet));
        }
    }
}
