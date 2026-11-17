//! Which mesh corners stand for their singular row (PK-1178; the
//! decision is the chart service's `corner_on_singular_row`).
//!
//! A triangle corner at a chart pole (a singular row's point, asked of
//! the chart service) is integrated as the quad from the pole
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
use popkern_core::geometry::chart::ChartDescriptor;
use popkern_core::geometry::chart::ops::{
    SINGULAR_EXACT, SINGULAR_RESOLVED, corner_on_singular_row, row_resolutions,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

/// A face's singular-row context: its loops' vertices, and per singular row of
/// its chart the distance within which a vertex is resolved at the row.
pub(super) struct FaceRows {
    vertices: Vec<FaceVertex>,
    resolved: Vec<f64>,
}

/// A vertex of a face's loops: its position and tolerance (its own,
/// never under the body's linear tolerance).
pub(super) struct FaceVertex {
    position: Point3,
    tolerance: f64,
}

/// `face`'s singular-row context over `descriptor` (its surface's chart).
pub(super) fn face_rows(
    brep: &BRep,
    face: FaceHandle,
    surface: Option<(&SurfaceGeometry, &ChartDescriptor)>,
) -> FaceRows {
    let resolved = surface.map_or_else(Vec::new, |(s, d)| row_resolutions(d, s, SINGULAR_RESOLVED));
    FaceRows {
        vertices: face_vertices(brep, face),
        resolved,
    }
}

/// The vertices of `face`'s loops.
fn face_vertices(brep: &BRep, face: FaceHandle) -> Vec<FaceVertex> {
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

/// The corner of a triangle at a chart pole or a cone's apex, and whether
/// its `u` (else its `v`) is the collapsed coordinate: asked of the chart
/// service in space, never by partial ratios. A corner stands for its
/// singular row when its point lies on the row's point to rounding (any
/// mesh point; a primitive's pole vertex is exact in space even where its
/// stored (u, v) sits a few micro-radians off), or when it is a face
/// vertex within its tolerance of the point (a boolean's apex vertex on a
/// kilometre cone sits within its tolerance of the apex, not at it,
/// PK-1133) and resolved at the row ([`SINGULAR_RESOLVED`] of the row's
/// distortion radius, or of the other corners' distance from it).
pub(super) fn singular_corner(
    descriptor: &ChartDescriptor,
    points: [Point3; 3],
    (rows, gate): (&FaceRows, f64),
) -> Option<(usize, bool)> {
    if descriptor.singular_rows.is_empty() {
        return None;
    }
    let scale = points
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(1.0, f64::max);
    (0..3).find_map(|k| {
        let tolerance = rows
            .vertices
            .iter()
            .find(|q| q.position.distance_to(&points[k]) <= gate)
            .map(|q| q.tolerance);
        let exact = SINGULAR_EXACT * scale;
        corner_on_singular_row(
            descriptor,
            &rows.resolved,
            points,
            k,
            (tolerance, exact, SINGULAR_RESOLVED),
        )
        .map(|(_, collapsed_u)| (k, collapsed_u))
    })
}
