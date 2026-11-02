//! Where along an edge its radial order is taken.
//!
//! The order is taken at a few stations and must agree at all of them:
//! two faces crossing along the edge (a tangency inside it) would swap
//! fins between stations, and the edge declines `Inconsistent` instead
//! of serving one station's guess. A curved edge gets five stations, a
//! straight one three. A straight edge whose every fin lies on a plane
//! gets one: each fin's inward direction n × t is then the same at
//! every point of the edge, so the order is too, and the other stations
//! only repeat it (on a 48-gon prism's union with a box, 0.24 of the
//! complex's 0.93 ms).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;

use super::ids::{EdgeId, FinId};
use super::model::CellComplex;
use super::radial::{FaceGeometry, station};

/// The stations (edge fractions) `edge`'s order is taken at, the first
/// deciding it. `length` is the edge's sampled length.
pub(super) fn stations(
    c: &CellComplex,
    geo: &FaceGeometry<'_>,
    edge: EdgeId,
    length: f64,
    linear: f64,
) -> &'static [f64] {
    let e = &c.edges[edge.index()];
    let chord = c.vertices[e.start.index()]
        .position
        .distance_to(&c.vertices[e.end.index()].position);
    if length > chord + linear {
        &[0.5, 0.25, 0.75, 0.125, 0.875]
    } else if e.radial.iter().all(|&fin| on_plane(c, geo, fin)) {
        &[0.5]
    } else {
        &[0.5, 0.25, 0.75]
    }
}

/// Whether `fin`'s face lies on a plane (its primary owner's surface).
fn on_plane(c: &CellComplex, geo: &FaceGeometry<'_>, fin: FinId) -> bool {
    c.fins
        .get(fin.index())
        .and_then(|f| c.faces.get(f.face.index()))
        .and_then(|face| face.owners.first())
        .and_then(|o| geo.operands.get(o.operand)?.faces.get(o.face))
        .and_then(|f| f.surface.as_ref())
        .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Plane { .. }))
}

/// Length of an edge, sampled (closed edges have no chord).
pub(super) fn edge_length(c: &CellComplex, edge: EdgeId) -> f64 {
    let pts: Vec<Point3> = (0..=8)
        .filter_map(|k| station(c, edge, f64::from(k) / 8.0).map(|(p, _)| p))
        .collect();
    pts.windows(2).map(|w| w[0].distance_to(&w[1])).sum()
}

/// Whether two cyclic orders are rotations of each other.
pub(super) fn same_cycle(a: &[(FinId, bool)], b: &[(FinId, bool)]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let Some(start) = b.iter().position(|x| x.0 == a[0].0) else {
        return false;
    };
    (0..a.len()).all(|k| a[k] == b[(start + k) % b.len()])
}
