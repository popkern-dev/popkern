//! Feature edges and blend roles from a solid's welded mesh and its
//! faces' exact surfaces.
//!
//! At every mesh edge where two faces meet, both faces' exact outward
//! normals are evaluated at the edge's two vertices (a freeform face's
//! through the caller's surface evaluator; its triangle's normal only
//! when that fails). The angle between them is the
//! edge's turning angle: below `tangent_angle` the faces join G1, above
//! `sharp_angle` the edge is sharp, and a sharp edge is concave when the
//! faces turn inward ((n_a × n_b) · e < 0 for the edge direction e as
//! traversed by face a's triangle). Concave sharp segments of one face
//! pair are chained into polylines.
//!
//! A cylinder or torus face that joins at least two other faces G1 is
//! inferred as a fillet of its (minor) radius. This also matches a full
//! round (the half-cylinder end of a slot); the role is `Inferred`. A
//! freeform face joined G1 to at least two faces (a variable-radius or
//! freeform blend) is a fillet too, with the smallest radius of
//! curvature its exact normals show along the tessellation's edges.
//!
//! A planar face that meets exactly two planar faces obliquely (each
//! sharp and under 85°, both convex or both concave), whose turning
//! angles add up to the angle between those two faces, and which is a
//! narrow strip, is a chamfer of the edge they would have met at; where
//! such candidates share an edge, the largest set of non-adjacent ones
//! are the chamfers.
//!
//! `popkern_measure::dihedral_angle` at the pin takes one normal per
//! face, which only suits planes; the per-vertex exact normals here
//! cover curved faces too.

use std::collections::BTreeMap;

use popkern_analysis::{ConcaveEdge, FaceId, FaceInfo, Surface};
use popkern_math::PortableMath;

use crate::geometry::{cross, dot, norm, sub, unit};

/// Angle thresholds for feature detection, in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FeatureTolerances {
    /// Largest normal angle at which two faces count as joining G1.
    pub tangent_angle: f64,
    /// Smallest normal angle at which an edge counts as sharp.
    pub sharp_angle: f64,
}

impl Default for FeatureTolerances {
    /// G1 within 1e-3 rad (an exact-geometry tangency measures about
    /// 1e-15 on the fixtures; imported CAD tangency is looser), sharp
    /// above 5°.
    fn default() -> Self {
        FeatureTolerances {
            tangent_angle: 1.0e-3,
            sharp_angle: 5.0 * std::f64::consts::PI / 180.0,
        }
    }
}

/// The exact outward-agnostic unit normal of a freeform face at a point
/// on it, when the surface can be evaluated there.
pub type NormalAt<'a> = dyn Fn(FaceId, [f64; 3]) -> Option<[f64; 3]> + 'a;

/// Features found on one solid.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SolidFeatures {
    /// Concave sharp edges, ordered by face pair then first vertex.
    pub concave_edges: Vec<ConcaveEdge>,
    /// Largest normal angle seen on any G1 face join (diagnostic for
    /// the tangent threshold).
    pub max_tangent_angle: f64,
}

/// Angle between unit vectors, accurate near 0 and π.
pub(crate) fn angle(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(cross(a, b)).patan2(dot(a, b))
}

/// The exact outward normal of `surface` at `p`, signed to agree with
/// the outward triangle normal `hint`; `None` for a freeform surface or
/// a point on an axis.
fn surface_normal(surface: &Surface, p: [f64; 3], hint: [f64; 3]) -> Option<[f64; 3]> {
    let radial = |origin: [f64; 3], axis: [f64; 3]| {
        let v = sub(p, origin);
        let a = dot(v, axis);
        (unit([0, 1, 2].map(|k| v[k] - a * axis[k])), a)
    };
    let n = match *surface {
        Surface::Plane { normal, .. } => Some(normal),
        Surface::Cylinder { origin, axis, .. } => radial(origin, axis).0,
        Surface::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let (r, a) = radial(apex, axis);
            let (s, c) = (half_angle.psin(), half_angle.pcos());
            let sign = if a < 0.0 { -1.0 } else { 1.0 };
            r.and_then(|r| unit([0, 1, 2].map(|k| r[k] * c - sign * s * axis[k])))
        }
        Surface::Sphere { centre, .. } => unit(sub(p, centre)),
        Surface::Torus {
            centre,
            axis,
            major_radius,
            ..
        } => radial(centre, axis).0.and_then(|r| {
            let tube = [0, 1, 2].map(|k| centre[k] + major_radius * r[k]);
            unit(sub(p, tube))
        }),
        Surface::Freeform => None,
    }?;
    Some(if dot(n, hint) < 0.0 { n.map(|c| -c) } else { n })
}

/// `n` turned to agree with `hint`.
fn signed(n: [f64; 3], hint: [f64; 3]) -> [f64; 3] {
    if dot(n, hint) < 0.0 { n.map(|c| -c) } else { n }
}

fn triangle_normal(v: &[[f64; 3]], t: [u32; 3]) -> Option<[f64; 3]> {
    let [a, b, c] = t.map(|i| v[i as usize]);
    unit(cross(sub(b, a), sub(c, a)))
}

/// One mesh edge between two faces: (face a, face b, a's edge start,
/// a's edge end, turning angle, concave).
struct Segment {
    faces: (FaceId, FaceId),
    ends: (u32, u32),
    turning: f64,
    concave: bool,
}

/// The face-boundary segments of a closed welded mesh.
fn boundary_segments(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    triangle_face: &[FaceId],
    surface_of: &dyn Fn(FaceId) -> Option<Surface>,
    freeform: &NormalAt<'_>,
) -> Vec<Segment> {
    // Directed edge (a, b) → triangle; its twin is (b, a).
    let mut directed: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for (i, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            directed.insert((t[k], t[(k + 1) % 3]), i);
        }
    }
    let mut out = Vec::new();
    for (&(a, b), &ta) in &directed {
        let Some(&tb) = directed.get(&(b, a)) else {
            continue;
        };
        let (fa, fb) = (triangle_face[ta], triangle_face[tb]);
        // Each undirected edge once, from the lower face's side.
        if fa >= fb {
            continue;
        }
        let (Some(ha), Some(hb)) = (
            triangle_normal(vertices, triangles[ta]),
            triangle_normal(vertices, triangles[tb]),
        ) else {
            continue;
        };
        let (sa, sb) = (surface_of(fa), surface_of(fb));
        let mut turning: f64 = 0.0;
        let mut turn = [0.0; 3];
        let exact = |face: FaceId, s: &Option<Surface>, p: [f64; 3], hint: [f64; 3]| match s {
            Some(Surface::Freeform) => freeform(face, p).map(|n| signed(n, hint)),
            Some(s) => surface_normal(s, p, hint),
            None => None,
        };
        for v in [a, b] {
            let p = vertices[v as usize];
            let na = exact(fa, &sa, p, ha).unwrap_or(ha);
            let nb = exact(fb, &sb, p, hb).unwrap_or(hb);
            turning = turning.max(angle(na, nb));
            let x = cross(na, nb);
            for k in 0..3 {
                turn[k] += x[k];
            }
        }
        let e = sub(vertices[b as usize], vertices[a as usize]);
        out.push(Segment {
            faces: (fa, fb),
            ends: (a, b),
            turning,
            concave: dot(turn, e) < 0.0,
        });
    }
    out
}

/// How one face pair joins: whether every shared segment is G1, its
/// concave sharp segments with their summed turning angle, and over all
/// segments their count, total length, summed turning and whether every
/// one is sharp and convex.
pub(crate) struct PairJoin {
    pub tangent: bool,
    pub concave: Vec<(u32, u32)>,
    pub turning_sum: f64,
    pub segments: usize,
    pub length: f64,
    pub all_turning: f64,
    pub sharp_convex: bool,
}

/// Chains segments (as vertex pairs) into maximal polylines, in a fixed
/// order: open chains start at their lowest free end.
fn chain(pairs: &[(u32, u32)]) -> Vec<Vec<u32>> {
    let mut adj: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (i, &(a, b)) in pairs.iter().enumerate() {
        adj.entry(a).or_default().push(i);
        adj.entry(b).or_default().push(i);
    }
    let mut used = vec![false; pairs.len()];
    let mut chains = Vec::new();
    // Free ends (degree 1) first, then whatever closed loops remain.
    let starts: Vec<u32> = adj
        .iter()
        .filter(|(_, s)| s.len() == 1)
        .map(|(&v, _)| v)
        .chain(pairs.iter().map(|p| p.0))
        .collect();
    for start in starts {
        let mut line = vec![start];
        let mut at = start;
        while let Some(&i) = adj.get(&at).and_then(|s| s.iter().find(|&&i| !used[i])) {
            used[i] = true;
            let (a, b) = pairs[i];
            at = if a == at { b } else { a };
            line.push(at);
        }
        if line.len() > 1 {
            chains.push(line);
        }
    }
    chains
}

/// Finds one solid's concave sharp edges and sets the inferred blend
/// roles (fillet, chamfer) of its faces. `faces` are this solid's faces
/// with consecutive ids; planar normals must already point outward.
/// `freeform` evaluates exact normals on freeform faces.
pub fn detect_features(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    triangle_face: &[FaceId],
    faces: &mut [FaceInfo],
    tolerances: &FeatureTolerances,
    freeform: &NormalAt<'_>,
) -> SolidFeatures {
    let first = faces.first().map_or(0, |f| f.id.0 as usize);
    let index = |f: FaceId| (f.0 as usize).checked_sub(first);
    let snapshot: Vec<Surface> = faces.iter().map(|f| f.surface.clone()).collect();
    let surface_of = |f: FaceId| index(f).and_then(|i| snapshot.get(i)).cloned();
    let segments = boundary_segments(vertices, triangles, triangle_face, &surface_of, freeform);

    let mut features = SolidFeatures::default();
    let mut pairs: BTreeMap<(FaceId, FaceId), PairJoin> = BTreeMap::new();
    for s in &segments {
        let entry = pairs.entry(s.faces).or_insert(PairJoin {
            tangent: true,
            concave: Vec::new(),
            turning_sum: 0.0,
            segments: 0,
            length: 0.0,
            all_turning: 0.0,
            sharp_convex: true,
        });
        entry.segments += 1;
        entry.length += norm(sub(
            vertices[s.ends.1 as usize],
            vertices[s.ends.0 as usize],
        ));
        entry.all_turning += s.turning;
        if s.turning < tolerances.sharp_angle || s.concave {
            entry.sharp_convex = false;
        }
        if s.turning > tolerances.tangent_angle {
            entry.tangent = false;
        } else {
            features.max_tangent_angle = features.max_tangent_angle.max(s.turning);
        }
        if s.turning >= tolerances.sharp_angle && s.concave {
            entry.concave.push(s.ends);
            entry.turning_sum += s.turning;
        }
    }
    for (&(fa, fb), join) in &pairs {
        if join.concave.is_empty() {
            continue;
        }
        let dihedral = join.turning_sum / join.concave.len() as f64;
        for line in chain(&join.concave) {
            features.concave_edges.push(ConcaveEdge {
                faces: (fa, fb),
                dihedral,
                polyline: line.iter().map(|&v| vertices[v as usize]).collect(),
            });
        }
    }
    crate::roles::infer_roles(
        &crate::roles::Mesh {
            vertices,
            triangles,
            triangle_face,
        },
        faces,
        &pairs,
        tolerances,
        freeform,
    );
    features
}
