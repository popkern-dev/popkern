//! A design's surface as a part (PH-117): an analysis model the
//! strength preview can solve, with the design's supports and loads.
//!
//! Each surface triangle lying on one of the original part's faces
//! (its centroid within `tolerance` of that face's triangles, normals
//! within ~37°) keeps that face's id, so the setup's named supports and
//! loads still find it; every other triangle (the design's new
//! surfaces) belongs to one more face, free. The new body is one solid.

use std::collections::HashMap;

use crate::model::{AnalysisModel, FaceId, FaceInfo, RoleSource, SolidId, SolidMesh, Surface};

use super::surface::SurfaceMesh;

/// Smallest cosine between a design triangle's normal and the original
/// face's for the triangle to belong to that face.
const NORMAL_COSINE: f64 = 0.8;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Unit normal and twice the area of triangle `t` of `vertices`.
fn normal(vertices: &[[f64; 3]], t: [u32; 3]) -> ([f64; 3], f64) {
    let [a, b, c] = t.map(|i| vertices[i as usize]);
    let n = cross(sub(b, a), sub(c, a));
    let len = dot(n, n).sqrt();
    if len > 0.0 {
        (n.map(|x| x / len), len)
    } else {
        ([0.0; 3], 0.0)
    }
}

/// The design body of `surface`, its faces inherited from `model`
/// within `tolerance` (mm). The new surfaces' face is the last,
/// `FaceId(model.faces.len())`.
pub fn design_model(model: &AnalysisModel, surface: &SurfaceMesh, tolerance: f64) -> AnalysisModel {
    let bucket = tolerance.max(f64::MIN_POSITIVE);
    let key = |p: [f64; 3]| p.map(|x| (x / bucket).floor() as i64);
    // The original triangles by bucket, in (solid, triangle) order;
    // looked up only, never iterated.
    let mut buckets: HashMap<[i64; 3], Vec<(usize, usize)>> = HashMap::new();
    for (s, solid) in model.solids.iter().enumerate() {
        for (t, tri) in solid.triangles.iter().enumerate() {
            let pts = tri.map(|i| solid.vertices[i as usize]);
            let lo = key([0, 1, 2]
                .map(|a| pts.iter().map(|p| p[a]).fold(f64::INFINITY, f64::min) - tolerance));
            let hi = key([0, 1, 2]
                .map(|a| pts.iter().map(|p| p[a]).fold(f64::NEG_INFINITY, f64::max) + tolerance));
            for z in lo[2]..=hi[2] {
                for y in lo[1]..=hi[1] {
                    for x in lo[0]..=hi[0] {
                        buckets.entry([x, y, z]).or_default().push((s, t));
                    }
                }
            }
        }
    }
    let design_face = FaceId(model.faces.len() as u32);
    let mut triangle_face = Vec::with_capacity(surface.triangles.len());
    for &t in &surface.triangles {
        let (n, area2) = normal(&surface.vertices, t);
        let pts = t.map(|i| surface.vertices[i as usize]);
        let c = [0, 1, 2].map(|a| (pts[0][a] + pts[1][a] + pts[2][a]) / 3.0);
        let mut best: Option<(f64, FaceId)> = None;
        if area2 > 0.0 {
            for &(s, o) in buckets.get(&key(c)).map_or(&[][..], Vec::as_slice) {
                let solid = &model.solids[s];
                let tri = solid.triangles[o];
                let (m, oa) = normal(&solid.vertices, tri);
                if oa <= 0.0 || dot(n, m) < NORMAL_COSINE {
                    continue;
                }
                let corners: Vec<[f64; 3]> =
                    tri.iter().map(|&i| solid.vertices[i as usize]).collect();
                let d = super::super::infill::polygon_distance(c, &corners, m);
                if d <= tolerance && best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, solid.triangle_face[o]));
                }
            }
        }
        triangle_face.push(best.map_or(design_face, |(_, f)| f));
    }
    let mut areas = vec![0.0; model.faces.len() + 1];
    for (&t, f) in surface.triangles.iter().zip(&triangle_face) {
        areas[f.0 as usize] += 0.5 * normal(&surface.vertices, t).1;
    }
    let mut faces: Vec<FaceInfo> = model
        .faces
        .iter()
        .map(|f| FaceInfo {
            solid: SolidId(0),
            surface: Surface::Freeform,
            role: None,
            role_source: RoleSource::Inferred,
            area: areas[f.id.0 as usize],
            exact_area: None,
            folded_area: None,
            ..f.clone()
        })
        .collect();
    faces.push(FaceInfo {
        id: design_face,
        solid: SolidId(0),
        name: Some("design surface".to_string()),
        surface: Surface::Freeform,
        role: None,
        role_source: RoleSource::Inferred,
        area: areas[design_face.0 as usize],
        exact_area: None,
        folded_area: None,
    });
    AnalysisModel {
        units: model.units,
        chord_tolerance: model.chord_tolerance,
        solids: vec![SolidMesh {
            id: SolidId(0),
            name: Some("design".to_string()),
            vertices: surface.vertices.clone(),
            triangles: surface.triangles.clone(),
            triangle_face,
            exact_volume: None,
        }],
        faces,
        concave_edges: Vec::new(),
        provenance: model.provenance.clone(),
    }
}
