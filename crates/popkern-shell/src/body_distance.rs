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
    /// Whether `polygons[0]` is the OUTER loop's polygon (false when
    /// the outer loop was skipped as degenerate).
    outer_polygon_first: bool,
    /// Parameter periodicity of the surface (query canonicalization).
    periodicity: (Option<f64>, Option<f64>),
}

/// Trim-polygon index over all faces of a body.
pub(crate) struct BodyDistanceIndex {
    pub entries: Vec<FaceTrimEntry>,
}

impl FaceTrimEntry {
    /// The OUTER loop's UV polygon in traversal order (`None` when the
    /// outer loop produced no usable polygon — inner-only entries
    /// cannot answer outer-winding questions).
    pub(crate) fn outer_polygon(&self) -> Option<&[(f64, f64)]> {
        if self.outer_polygon_first {
            self.polygons.first().map(|p| p.as_slice())
        } else {
            None
        }
    }
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

/// Nearest point of a trimmed face patch to a query point.
pub(crate) struct PatchNearest {
    /// The nearest point on (or near, for the sampled-boundary arm)
    /// the trimmed patch.
    pub point: Point3,
    /// UV of the nearest point: exact for the surface arm, the nearest
    /// boundary SAMPLE's uv for the boundary arm (a projection seed,
    /// not an exact preimage).
    pub uv: (f64, f64),
    /// Distance from the query point.
    pub distance: f64,
}

impl FaceTrimEntry {
    /// Build the trim entry for a single face (`None` when the face
    /// has no surface or no walkable loop with ≥3 samples).
    pub(crate) fn build(body: &BRep, fh: FaceHandle) -> Option<FaceTrimEntry> {
        build_entry(body, fh)
    }

    /// Distance from `p` to this trimmed face patch.
    pub(crate) fn distance(&self, p: &Point3) -> Option<f64> {
        self.nearest(p).map(|n| n.distance)
    }

    /// Nearest point of this trimmed face patch to `p` (PK-328
    /// clearance queries); same arms and comparison order as the
    /// original distance query.
    pub(crate) fn nearest(&self, p: &Point3) -> Option<PatchNearest> {
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
        let mut boundary_point = self.boundary[nearest_idx];
        for j in [
            nearest_idx.saturating_sub(1),
            nearest_idx,
            (nearest_idx + 1).min(self.boundary.len() - 1),
        ] {
            if j + 1 < self.boundary.len() {
                let (d, q) = segment_nearest(p, &self.boundary[j], &self.boundary[j + 1]);
                if d < boundary_arm {
                    boundary_arm = d;
                    boundary_point = q;
                }
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
            let d = p.distance_to(&q);
            if d < boundary_arm {
                return Some(PatchNearest {
                    point: q,
                    uv: (u, v),
                    distance: d,
                });
            }
        }
        Some(PatchNearest {
            point: boundary_point,
            uv: seed,
            distance: boundary_arm,
        })
    }

    /// The sampled boundary points (trim-true clearance sweeps,
    /// PK-328).
    pub(crate) fn boundary_samples(&self) -> &[Point3] {
        &self.boundary
    }

    /// Unit surface normal at a boundary sample index (`None` at
    /// degenerate points).
    pub(crate) fn boundary_normal(&self, idx: usize) -> Option<popkern_core::geometry::Vec3> {
        let uv = *self.boundary_uv.get(idx)?;
        self.surface_normal(uv)
    }

    /// Unit surface normal at a UV (`None` at degenerate points).
    pub(crate) fn surface_normal(&self, uv: (f64, f64)) -> Option<popkern_core::geometry::Vec3> {
        let n = self.surface.normal_at(uv.0, uv.1);
        let len = n.length();
        (len > 1e-12).then(|| n * (1.0 / len))
    }

    /// Even-odd containment of a UV point over all loop polygons.
    pub(crate) fn contains(&self, u: f64, v: f64) -> bool {
        // Canonicalize the query onto each polygon's period branch:
        // polygons are chain-unwrapped at build time (a seam-straddling
        // sphere octant's raw pcurves carry mixed branches — even-odd
        // over the broken shape accepted the WHOLE sphere and a
        // preview leaked out-of-trim samples 0.09 deep, PK-328), so a
        // raw-branch query must be shifted to match.
        let mut inside = false;
        for polygon in &self.polygons {
            if polygon.is_empty() {
                continue;
            }
            // Anchor at the polygon's bbox CENTER: a full-wrap ring
            // spans a whole period, and canonicalizing to its first
            // point ejects queries near the far end of the span.
            let (mut ulo, mut uhi, mut vlo, mut vhi) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
            for &(pu, pv) in polygon {
                ulo = ulo.min(pu);
                uhi = uhi.max(pu);
                vlo = vlo.min(pv);
                vhi = vhi.max(pv);
            }
            let uq = canonicalize(u, (ulo + uhi) * 0.5, self.periodicity.0);
            let vq = canonicalize(v, (vlo + vhi) * 0.5, self.periodicity.1);
            if polygon_contains(polygon, uq, vq) {
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
    let periodicity = surface.parameter_periodicity();
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());

    let mut polygons = Vec::new();
    let mut boundary = Vec::new();
    let mut boundary_uv = Vec::new();
    let mut outer_polygon_first = false;
    for (loop_idx, lh) in loops.into_iter().enumerate() {
        let lh: popkern_core::topology::LoopHandle = lh;
        let is_outer = loop_idx == 0;
        let loop_data = body.loops.get(lh)?;
        let start = loop_data.half_edge;
        if !body.half_edges.contains_key(start) {
            continue;
        }
        let mut polygon: Vec<(f64, f64)> = Vec::new();
        // Chain-unwrapped push: successive samples stay on one period
        // branch (see `FaceTrimEntry::contains`).
        let push_unwrapped = |polygon: &mut Vec<(f64, f64)>, u: f64, v: f64| {
            let (u, v) = match polygon.last() {
                Some(&(pu, pv)) => (
                    canonicalize(u, pu, periodicity.0),
                    canonicalize(v, pv, periodicity.1),
                ),
                None => (u, v),
            };
            polygon.push((u, v));
        };
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
                    push_unwrapped(&mut polygon, q.x, q.y);
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
                    push_unwrapped(&mut polygon, uv.0, uv.1);
                }
            }
            cursor = he.next;
            steps += 1;
            if cursor == start || steps > bound {
                break;
            }
        }
        if polygon.len() >= 3 {
            if is_outer {
                outer_polygon_first = true;
            }
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
        outer_polygon_first,
        periodicity,
    })
}

/// Shift `x` by whole periods to land within half a period of `anchor`.
fn canonicalize(x: f64, anchor: f64, period: Option<f64>) -> f64 {
    match period {
        Some(p) if p > 0.0 => {
            let mut x = x;
            while x - anchor > p * 0.5 {
                x -= p;
            }
            while anchor - x > p * 0.5 {
                x += p;
            }
            x
        }
        _ => x,
    }
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

fn segment_nearest(p: &Point3, a: &Point3, b: &Point3) -> (f64, Point3) {
    let ab = *b - *a;
    let len_sq = ab.dot(&ab);
    if len_sq <= 0.0 {
        return (p.distance_to(a), *a);
    }
    let t = ((*p - *a).dot(&ab) / len_sq).clamp(0.0, 1.0);
    let q = *a + ab * t;
    (p.distance_to(&q), q)
}
