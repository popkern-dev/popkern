//! Which pole corners stand for their pole row (PK-1178).
//!
//! A triangle corner at a chart pole (`pole_corner`: one partial
//! vanishing beside the other) is integrated as the quad from the pole
//! row under the other two corners (`surface_term`): a loop meeting the
//! pole meets it at a vertex, which the mesh keeps as one (u, v) point
//! whose collapsed coordinate is arbitrary, so its (u, v) triangles
//! cover one row of the fan their facets span. That holds for a vertex
//! at the pole within its tolerance (a primitive's pole, whose stored
//! (u, v) can sit 3e-6 rad off it; a boolean's apex vertex 0.016 up a
//! kilometre cone, PK-1133). A trim sample or a vertex beside a pole
//! carries its own longitude: the mesher tiles the chart with it
//! (`chart_degenerate` keeps it regular), its (u, v) triangles are
//! their patches, and the quad from the pole row under its neighbours
//! covers a sector their triangles cover too. A ball cut by a plane
//! whose circle encircles the pole 1e-5 from it, turned, read its ring
//! samples and its seam crossing there as poles and its volume 4.5e-5
//! high (one triangle 8.0e-4 against a facet of 2.2e-6).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle};

/// A collapsed partial under this share of the other is the pole to
/// rounding, whatever lies there.
const POLE_EXACT: f64 = 1e-12;

/// A vertex of a face's loops: its position and tolerance (its own,
/// never under the body's linear tolerance).
pub(super) struct FaceVertex {
    position: Point3,
    tolerance: f64,
}

/// The vertices of `face`'s loops.
pub(super) fn face_vertices(brep: &BRep, face: FaceHandle) -> Vec<FaceVertex> {
    let Some(data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let linear = brep.tolerance.linear();
    let mut out = Vec::new();
    for &lp in std::iter::once(&data.outer_loop).chain(&data.inner_loops) {
        let Some(start) = brep.loops.get(lp).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        while let Some(h) = brep.half_edges.get(he) {
            if let Some(v) = brep.vertices.get(h.origin) {
                out.push(FaceVertex {
                    position: v.position,
                    tolerance: v.tolerance.unwrap_or(0.0).max(linear),
                });
            }
            he = h.next;
            if he == start || out.len() > 100_000 {
                break;
            }
        }
    }
    out
}

/// The pole corner of `uv` (`super::pole_corner`) that stands for its
/// pole row: a corner whose collapsed partial vanishes to rounding, or
/// one at a face vertex (within `gate`) lying at the pole within the
/// vertex's tolerance, its distance off the axis read from its position
/// (half the move of half a period of the collapsed coordinate there).
/// Every corner read as a pole is asked: a ring sample beside a pole
/// and a pole-row sample in one triangle stand for the row by the
/// latter. A collapsed coordinate that is not periodic keeps the
/// ratio's reading.
pub(super) fn mesh_pole(
    surface: &SurfaceGeometry,
    uv: [(f64, f64); 3],
    points: [Point3; 3],
    (vertices, gate): (&[FaceVertex], f64),
) -> Option<(usize, bool)> {
    let (pu, pv) = surface.parameter_periodicity();
    let row = |k: usize, collapsed_u: bool| {
        let (u, v) = uv[k];
        let (su, sv) = (
            surface.partial_u(u, v).length(),
            surface.partial_v(u, v).length(),
        );
        if (collapsed_u && su <= POLE_EXACT * sv) || (!collapsed_u && sv <= POLE_EXACT * su) {
            return true;
        }
        let Some(vertex) = vertices
            .iter()
            .find(|q| q.position.distance_to(&points[k]) <= gate)
        else {
            return false;
        };
        let period = if collapsed_u { pu } else { pv };
        let Some(p) = period.filter(|&p| p > 0.0) else {
            return true;
        };
        let Some((u, v)) = surface.closest_point(&points[k]) else {
            return false;
        };
        let shifted = if collapsed_u {
            surface.point_at(u + 0.5 * p, v)
        } else {
            surface.point_at(u, v + 0.5 * p)
        };
        0.5 * shifted.distance_to(&points[k]) <= vertex.tolerance
    };
    super::pole_corner(surface, uv, &row)
}
