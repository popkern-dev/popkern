//! The window's boundary on the branch its edges lie on in space (W2c,
//! chart-service.md §6b step 2).
//!
//! A boundary pcurve on a periodic chart can be stored the long way round
//! its axis: its ends agree with its edge's, its interior runs the other
//! way about (a strip through the seam whose bottom edge was written u 5.5
//! → 0.8 rather than 5.5 → 7.08). Lifting the record by whole periods
//! cannot meet such a pcurve where the edge is; the production finder's
//! half-period shift healed it by tearing the boundary. The window reads
//! the pcurve against its edge in space instead: where its interior lies
//! far from the edge while its ends lie on it, the boundary is the edge's
//! curve projected onto the surface and unwrapped continuously from the
//! pcurve's start (the chart service's nearest period offset), on the
//! pcurve's own parameter range.

use popkern_core::geometry::chart::ChartDescriptor;
use popkern_core::geometry::chart::ops::nearest_offset;
use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval as _, SurfaceGeometry};
use popkern_core::geometry::{Point3, chart::Axis};
use popkern_core::topology::{BRep, HalfEdgeHandle};

/// Interior fractions at which a pcurve is read against its edge.
const PROBES: [f64; 3] = [0.25, 0.5, 0.75];
/// Samples of an edge's curve for a re-lifted boundary.
const RELIFT_SAMPLES: usize = 65;

/// `boundary` with every pcurve stored the long way round re-lifted from
/// its edge (see the module docs); a pcurve with no edge curve, or whose
/// ends do not lie on its edge, is kept as stored.
pub(super) fn on_its_edges(
    brep: &BRep,
    surface: &SurfaceGeometry,
    boundary: &[(HalfEdgeHandle, NurbsCurve)],
) -> Vec<(HalfEdgeHandle, NurbsCurve)> {
    let descriptor = ChartDescriptor::of(surface, &brep.tolerance);
    boundary
        .iter()
        .map(|(h, pc)| {
            let lifted = relift(brep, (surface, &descriptor), *h, pc);
            (*h, lifted.unwrap_or_else(|| pc.clone()))
        })
        .collect()
}

/// The edge's curve point at fraction `f` of the half-edge's walk.
fn edge_point(brep: &BRep, h: HalfEdgeHandle) -> Option<impl Fn(f64) -> Point3 + '_> {
    let he = brep.half_edges.get(h)?;
    let edge = brep.edges.get(he.edge)?;
    let curve = edge.curve.as_ref()?;
    let (t0, t1) = if edge.half_edge == h {
        (curve.t_start, curve.t_end)
    } else {
        (curve.t_end, curve.t_start)
    };
    (t0.is_finite() && t1.is_finite())
        .then_some(move |f: f64| curve.geometry.point_at(t0 + (t1 - t0) * f))
}

fn relift(
    brep: &BRep,
    (surface, descriptor): (&SurfaceGeometry, &ChartDescriptor),
    h: HalfEdgeHandle,
    pc: &NurbsCurve,
) -> Option<NurbsCurve> {
    let period = descriptor.period(Axis::U).map(|p| p.period)?;
    let on_edge = edge_point(brep, h)?;
    let (a, b) = pc.parameter_range();
    if !(a.is_finite() && b.is_finite()) || b <= a {
        return None;
    }
    let on_surface = |f: f64| {
        let q = pc.point_at(a + (b - a) * f);
        surface.point_at(q.x, q.y)
    };
    // The edge's extent in space sets the scale: its ends agree within a
    // hundredth of it, and an interior a tenth of it away runs the other
    // way about (the long and short ways part by the diameter).
    let extent = (0..=8)
        .map(|k| on_edge(k as f64 / 8.0))
        .collect::<Vec<_>>()
        .windows(2)
        .map(|w| w[0].distance_to(&w[1]))
        .sum::<f64>()
        .max(brep.tolerance.linear());
    let ends_agree = [0.0, 1.0]
        .iter()
        .all(|&f| on_surface(f).distance_to(&on_edge(f)) <= 1e-2 * extent);
    let interior_apart = PROBES
        .iter()
        .any(|&f| on_surface(f).distance_to(&on_edge(f)) > 0.1 * extent);
    if !(ends_agree && interior_apart) {
        return None;
    }
    // The edge projected and unwrapped from the pcurve's start.
    let start = pc.point_at(a);
    let mut prev = (start.x, start.y);
    let mut points = Vec::with_capacity(RELIFT_SAMPLES);
    for k in 0..RELIFT_SAMPLES {
        let f = k as f64 / (RELIFT_SAMPLES - 1) as f64;
        let (u, v) = surface.closest_point(&on_edge(f))?;
        let u = prev.0 + nearest_offset(u, prev.0, period);
        prev = (u, v);
        points.push(Point3::new(u, v, 0.0));
    }
    let n = points.len();
    let mut knots = vec![a, a];
    knots.extend((1..n - 1).map(|k| a + (b - a) * k as f64 / (n - 1) as f64));
    knots.extend([b, b]);
    Some(NurbsCurve {
        degree: 1,
        knots,
        control_points: points,
        weights: vec![1.0; n],
    })
}
