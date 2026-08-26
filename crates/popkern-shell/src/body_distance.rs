//! Trimmed distance queries against a validated body's boundary
//! (PK-326 face-collapse detection; the same query underpins
//! `max_safe_offset` clearance work, PK-328).
//!
//! `distance` measures min over faces of the distance to the TRIMMED
//! face patch: the exact surface projection when its UV image lies
//! inside the face's trim polygons (built from the half-edge pcurves,
//! even-odd over outer + inner loops), else the distance to the
//! sampled boundary. Untrimmed surface distance would under-estimate
//! (an infinite plane passes arbitrarily close to points its face
//! never reaches) and silently kill live faces.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, uv_projection};
use popkern_core::topology::{BRep, FaceHandle};

/// Per-face trim data for distance queries.
pub(crate) struct FaceTrimEntry {
    /// The face this entry describes.
    pub face: FaceHandle,
    /// Face surface geometry.
    pub surface: SurfaceGeometry,
    /// Closed UV polygons, one per loop (outer first), chart-coherent
    /// per polygon.
    pub polygons: Vec<Vec<(f64, f64)>>,
    /// 3D images of every polygon sample (boundary samples).
    boundary: Vec<Point3>,
    /// UV of each boundary sample (projection branch seeds).
    boundary_uv: Vec<(f64, f64)>,
}

/// Trim-polygon index over all faces of a body.
pub(crate) struct BodyDistanceIndex {
    pub entries: Vec<FaceTrimEntry>,
}

const PCURVE_SAMPLES: usize = 8;

impl BodyDistanceIndex {
    /// Build the index. Faces without surfaces or walkable loops are
    /// skipped (they cannot bound distance).
    pub(crate) fn build(body: &BRep) -> Self {
        let mut entries = Vec::new();
        let mut face_handles: Vec<FaceHandle> = body.faces.keys().collect();
        face_handles.sort();
        for fh in face_handles {
            if let Some(entry) = build_entry(body, fh) {
                entries.push(entry);
            }
        }
        Self { entries }
    }

    /// Distance from `p` to the body's boundary (min over trimmed
    /// faces). `None` when no face produced a measurable arm.
    pub(crate) fn distance(&self, p: &Point3) -> Option<f64> {
        let mut best: Option<f64> = None;
        for entry in &self.entries {
            let Some(d) = entry.distance(p) else { continue };
            best = Some(best.map_or(d, |b: f64| b.min(d)));
        }
        best
    }
}

impl FaceTrimEntry {
    /// Build the trim entry for a single face (`None` when the face
    /// has no surface or no walkable loop with ≥3 samples).
    pub(crate) fn build(body: &BRep, fh: FaceHandle) -> Option<FaceTrimEntry> {
        build_entry(body, fh)
    }

    /// Distance from `p` to this trimmed face patch.
    fn distance(&self, p: &Point3) -> Option<f64> {
        // Boundary arm: nearest sample, refined by projection onto the
        // adjacent sample chords.
        let mut boundary_arm = f64::INFINITY;
        let mut nearest_idx = 0usize;
        for (i, b) in self.boundary.iter().enumerate() {
            let d = p.distance_to(b);
            if d < boundary_arm {
                boundary_arm = d;
                nearest_idx = i;
            }
        }
        if self.boundary.is_empty() {
            return None;
        }
        for j in [
            nearest_idx.saturating_sub(1),
            nearest_idx,
            (nearest_idx + 1).min(self.boundary.len() - 1),
        ] {
            if j + 1 < self.boundary.len() {
                boundary_arm = boundary_arm.min(segment_distance(
                    p,
                    &self.boundary[j],
                    &self.boundary[j + 1],
                ));
            }
        }

        // Surface arm: exact projection, accepted only inside the trim.
        let seed = self.boundary_uv[nearest_idx];
        let uv = uv_projection::project_point_to_uv_near(&self.surface, p, seed.0, seed.1)
            .or_else(|| uv_projection::project_point_to_uv(&self.surface, p));
        if let Some((u, v)) = uv
            && self.contains(u, v)
        {
            let q = self.surface.point_at(u, v);
            return Some(boundary_arm.min(p.distance_to(&q)));
        }
        Some(boundary_arm)
    }

    /// Even-odd containment of a UV point over all loop polygons.
    pub(crate) fn contains(&self, u: f64, v: f64) -> bool {
        let mut inside = false;
        for polygon in &self.polygons {
            if polygon_contains(polygon, u, v) {
                inside = !inside;
            }
        }
        inside
    }

    /// UV bounding box over every polygon.
    pub(crate) fn uv_bbox(&self) -> Option<(f64, f64, f64, f64)> {
        let mut it = self.polygons.iter().flatten();
        let first = it.next()?;
        let (mut u0, mut v0, mut u1, mut v1) = (first.0, first.1, first.0, first.1);
        for (u, v) in self.polygons.iter().flatten() {
            u0 = u0.min(*u);
            v0 = v0.min(*v);
            u1 = u1.max(*u);
            v1 = v1.max(*v);
        }
        Some((u0, v0, u1, v1))
    }
}

fn build_entry(body: &BRep, fh: FaceHandle) -> Option<FaceTrimEntry> {
    let face = body.faces.get(fh)?;
    let surface = face.surface.as_ref()?.geometry.clone();
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());

    let mut polygons = Vec::new();
    let mut boundary = Vec::new();
    let mut boundary_uv = Vec::new();
    for lh in loops {
        let loop_data = body.loops.get(lh)?;
        let start = loop_data.half_edge;
        if !body.half_edges.contains_key(start) {
            continue;
        }
        let mut polygon: Vec<(f64, f64)> = Vec::new();
        let mut cursor = start;
        let bound = body.half_edges.len();
        let mut steps = 0usize;
        loop {
            let he = body.half_edges.get(cursor)?;
            if let Some(pcurve) = he.pcurve.as_ref() {
                let (t0, t1) = pcurve.parameter_range();
                for i in 0..PCURVE_SAMPLES {
                    let t = t0 + (t1 - t0) * i as f64 / PCURVE_SAMPLES as f64;
                    let q = pcurve.point_at(t);
                    polygon.push((q.x, q.y));
                }
            } else {
                // Pcurve-less half-edge: project the 3D curve span,
                // branch-seeded from the previous polygon point.
                let edge = body.edges.get(he.edge)?;
                let curve = edge.curve.as_ref()?;
                use popkern_core::geometry::curve::CurveEval;
                for i in 0..PCURVE_SAMPLES {
                    let t = curve.t_start
                        + (curve.t_end - curve.t_start) * i as f64 / PCURVE_SAMPLES as f64;
                    let q3 = curve.geometry.point_at(t);
                    let uv = match polygon.last() {
                        Some(&(pu, pv)) => {
                            uv_projection::project_point_to_uv_near(&surface, &q3, pu, pv)
                        }
                        None => uv_projection::project_point_to_uv(&surface, &q3),
                    }?;
                    polygon.push(uv);
                }
            }
            cursor = he.next;
            steps += 1;
            if cursor == start || steps > bound {
                break;
            }
        }
        if polygon.len() >= 3 {
            for &(u, v) in &polygon {
                boundary.push(surface.point_at(u, v));
                boundary_uv.push((u, v));
            }
            polygons.push(polygon);
        }
    }
    if polygons.is_empty() {
        return None;
    }
    Some(FaceTrimEntry {
        face: fh,
        surface,
        polygons,
        boundary,
        boundary_uv,
    })
}

/// Standard even-odd ray cast (pnpoly).
fn polygon_contains(polygon: &[(f64, f64)], u: f64, v: f64) -> bool {
    let mut inside = false;
    let n = polygon.len();
    let mut j = n - 1;
    for i in 0..n {
        let (ui, vi) = polygon[i];
        let (uj, vj) = polygon[j];
        if ((vi > v) != (vj > v)) && (u < (uj - ui) * (v - vi) / (vj - vi) + ui) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn segment_distance(p: &Point3, a: &Point3, b: &Point3) -> f64 {
    let ab = *b - *a;
    let len_sq = ab.dot(&ab);
    if len_sq <= 0.0 {
        return p.distance_to(a);
    }
    let t = ((*p - *a).dot(&ab) / len_sq).clamp(0.0, 1.0);
    let q = *a + ab * t;
    p.distance_to(&q)
}
