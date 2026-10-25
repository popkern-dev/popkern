//! Tessellation folds: triangles turned against their face's surface.
//!
//! The welded mesh is watertight, so each face's tessellation is
//! bounded by the face's true edges; it can still fold over itself.
//! Each triangle's normal is compared with popkern's exact surface
//! normal at its centroid (projected to the surface), and the area the
//! face's minority sense covers is the folded area. A projection that
//! lands more than `max_distance` from the centroid found another part
//! of the surface (a helix's next turn: 72% of coil_spring's triangles)
//! and is not counted.

use std::collections::BTreeMap;

use popkern_analysis::FaceId;
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;

/// A triangle agrees with the surface when the cosine between its
/// normal and the exact normal is above this, and opposes it below the
/// negative; between, it is not counted.
const SENSE_COS: f64 = 0.5;

/// Folded area per face of `geometry` (the minority-sense area of its
/// triangles whose projection lands within `max_distance`).
pub(crate) fn folded_areas(
    geometry: &BTreeMap<FaceId, SurfaceGeometry>,
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    triangle_face: &[FaceId],
    max_distance: f64,
) -> BTreeMap<FaceId, f64> {
    // (area agreeing, area opposing) with the surface's own normal.
    let mut sense: BTreeMap<FaceId, (f64, f64)> = BTreeMap::new();
    for (t, f) in triangles.iter().zip(triangle_face) {
        let Some(g) = geometry.get(f) else {
            continue;
        };
        let entry = sense.entry(*f).or_insert((0.0, 0.0));
        let [a, b, c] = t.map(|i| vertices[i as usize]);
        let (u, v) = (
            [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
            [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
        );
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let twice = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if twice == 0.0 {
            continue;
        }
        let centroid = [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0);
        let p = Point3::new(centroid[0], centroid[1], centroid[2]);
        let Some((pu, pv)) = project_point_to_uv(g, &p) else {
            continue;
        };
        let q = g.point_at(pu, pv);
        let d = [q.x - p.x, q.y - p.y, q.z - p.z];
        if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() > max_distance {
            continue;
        }
        let e = g.normal_at(pu, pv);
        let len = (e.x * e.x + e.y * e.y + e.z * e.z).sqrt();
        if len == 0.0 {
            continue;
        }
        let cos = (n[0] * e.x + n[1] * e.y + n[2] * e.z) / (twice * len);
        if cos > SENSE_COS {
            entry.0 += 0.5 * twice;
        } else if cos < -SENSE_COS {
            entry.1 += 0.5 * twice;
        }
    }
    sense
        .into_iter()
        .map(|(f, (agree, oppose))| (f, agree.min(oppose)))
        .collect()
}

#[cfg(test)]
#[path = "fold_check_tests.rs"]
mod tests;
