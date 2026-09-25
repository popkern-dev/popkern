//! The setback boundary's sides (PK-527): one per strip, plus a bridge
//! on each reflex face corner of a saddle (see `junctions`). This
//! covers the side paths the extractor marches, the bridge's
//! cross-curvature from the neighbouring strips, and the periodic
//! branch rule for strip section paths.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{Point3, SurfaceGeometry, Vec3};

use super::junctions::{BridgeGeometry, JunctionKind, bridge_between, classify_junction};
use crate::fillet_surface::FilletSurfaceResult;

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
    /// A rail-flat bridge on its support plane.
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
