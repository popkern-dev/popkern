//! Fidelity gate for boundary-first meshes of geometrically closed
//! NURBS charts (extracted from `boundary_first.rs` under the
//! 500-line file budget). See [`closed_chart_fidelity_ok`].

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::FaceHandle;

use crate::TessellationOptions;
use crate::mesh::MeshVertex;

/// Whether a boundary-first result on a geometrically closed NURBS
/// chart is a faithful picture of the surface. A boundary-only
/// Delaunay is legal in UV yet can span the solid with 3D chords — a
/// closed helical tube's trim rails are its two seam chains, a full
/// profile period apart, so every triangle between them was a
/// tube-diameter chord and the coil meshed as a zero-volume ribbon
/// (PK-441). The MEDIAN triangle-edge midpoint deviation from the
/// surface gates at 4 × chord_height (the median separates a globally
/// degenerate sheet from a locally coarse spot); failures decline to
/// the uniform-grid path, recorded under the boundary-first trace
/// topic. Scoped to closed NURBS charts: several fillet-corpus faces
/// carry historically accepted coarse boundary-first meshes whose
/// neighbours conform to them — the general fidelity ratchet is the
/// chart-aware meshing campaign's scope (PK-445).
pub(super) fn closed_chart_fidelity_ok(
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    vertices: &[MeshVertex],
    triangles: &[[usize; 3]],
    options: &TessellationOptions,
) -> bool {
    let closed_chart = {
        let (pu, pv) = popkern_core::geometry::chart_periods::geometric_chart_periods(surface);
        matches!(surface, SurfaceGeometry::Nurbs(_)) && (pu.is_some() || pv.is_some())
    };
    if !closed_chart {
        return true;
    }
    let mut devs: Vec<f64> = Vec::with_capacity(triangles.len() * 3);
    for t in triangles {
        for (i, j) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let (a, b) = (&vertices[i], &vertices[j]);
            let um = 0.5 * (a.uv.0 + b.uv.0);
            let vm = 0.5 * (a.uv.1 + b.uv.1);
            let on_surf = surface.point_at(um, vm);
            let chord_mid = Point3::new(
                0.5 * (a.position.x + b.position.x),
                0.5 * (a.position.y + b.position.y),
                0.5 * (a.position.z + b.position.z),
            );
            devs.push(on_surf.distance_to(&chord_mid));
        }
    }
    devs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let median = devs.get(devs.len() / 2).copied().unwrap_or(0.0);
    if median > options.chord_height * 4.0 {
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!(
                "    [bf] face {face_handle:?} fidelity decline: median edge-mid deviation \
                 {median:.3e} > 4×chord_height"
            );
        }
        return false;
    }
    true
}

/// The `[bf-vol]` diagnostic under the boundary-first-ring topic:
/// wedge volume, UV area, and ring area of a produced face mesh.
pub(super) fn trace_bf_volume(
    face_handle: FaceHandle,
    vertices: &[MeshVertex],
    triangles: &[[usize; 3]],
    rings: &[super::BoundaryRing],
) {
    if !popkern_core::trace::trace_enabled("boundary-first-ring") {
        return;
    }
    let mut vol6 = 0.0;
    for t in triangles {
        let (a, b, c) = (
            vertices[t[0]].position,
            vertices[t[1]].position,
            vertices[t[2]].position,
        );
        vol6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    let mut uv_area2 = 0.0;
    let mut ring_area2 = 0.0_f64;
    for t in triangles {
        let (a, b, c) = (vertices[t[0]].uv, vertices[t[1]].uv, vertices[t[2]].uv);
        uv_area2 += ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
    }
    for ring in rings {
        let m = ring.uv.len();
        let mut s = 0.0;
        for i in 0..m {
            let (x1, y1) = ring.uv[i];
            let (x2, y2) = ring.uv[(i + 1) % m];
            s += x1 * y2 - x2 * y1;
        }
        ring_area2 += s;
    }
    eprintln!(
        "    [bf-vol] face {face_handle:?} tris {} wedge6 {vol6:.1} uv_area2 {uv_area2:.2} ring_area2 {:.2}",
        triangles.len(),
        ring_area2.abs()
    );
}
