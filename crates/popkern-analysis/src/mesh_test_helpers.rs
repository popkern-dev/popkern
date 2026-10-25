//! Closed, outward test meshes with known geometry.

use popkern_math::PortableMath;

use crate::model::{FaceId, SolidId, SolidMesh};

/// An axis-aligned box [min, max], one face id per side (0: −x, 1: +x,
/// 2: −y, 3: +y, 4: −z, 5: +z), two triangles per side.
pub(crate) fn box_mesh(min: [f64; 3], max: [f64; 3]) -> SolidMesh {
    let corner = |i: usize| {
        [
            if i & 1 == 0 { min[0] } else { max[0] },
            if i & 2 == 0 { min[1] } else { max[1] },
            if i & 4 == 0 { min[2] } else { max[2] },
        ]
    };
    let vertices: Vec<[f64; 3]> = (0..8).map(corner).collect();
    // Quads wound counter-clockwise seen from outside.
    let quads: [[u32; 4]; 6] = [
        [0, 4, 6, 2],
        [1, 3, 7, 5],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 2, 3, 1],
        [4, 5, 7, 6],
    ];
    let mut triangles = Vec::new();
    let mut triangle_face = Vec::new();
    for (f, q) in quads.iter().enumerate() {
        triangles.push([q[0], q[1], q[2]]);
        triangles.push([q[0], q[2], q[3]]);
        triangle_face.extend([FaceId(f as u32); 2]);
    }
    SolidMesh {
        id: SolidId(0),
        name: None,
        vertices,
        triangles,
        triangle_face,
        exact_volume: Some((0..3).map(|a| max[a] - min[a]).product()),
    }
}

/// Rotates every vertex by `degrees` about the axis (1, 2, 3)/√14
/// through `centre` (Rodrigues).
pub(crate) fn rotated(mut mesh: SolidMesh, degrees: f64, centre: [f64; 3]) -> SolidMesh {
    let n = 14.0f64.sqrt();
    let k = [1.0 / n, 2.0 / n, 3.0 / n];
    let theta = degrees * std::f64::consts::PI / 180.0;
    let (s, c) = (theta.psin(), theta.pcos());
    for v in &mut mesh.vertices {
        let p = [v[0] - centre[0], v[1] - centre[1], v[2] - centre[2]];
        let kxp = [
            k[1] * p[2] - k[2] * p[1],
            k[2] * p[0] - k[0] * p[2],
            k[0] * p[1] - k[1] * p[0],
        ];
        let kp = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
        for a in 0..3 {
            v[a] = centre[a] + p[a] * c + kxp[a] * s + k[a] * kp * (1.0 - c);
        }
    }
    mesh
}

/// A geodesic sphere: an icosahedron subdivided `levels` times with
/// vertices projected to radius `r` about `centre`; one face id.
pub(crate) fn icosphere(centre: [f64; 3], r: f64, levels: usize) -> SolidMesh {
    let t = (1.0 + 5.0f64.sqrt()) / 2.0;
    let mut v: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut f: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    let unit = |p: [f64; 3]| {
        let l = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        p.map(|c| c / l)
    };
    v = v.into_iter().map(unit).collect();
    for _ in 0..levels {
        let mut mids: Vec<((u32, u32), u32)> = Vec::new();
        let mut mid = |a: u32, b: u32, v: &mut Vec<[f64; 3]>| -> u32 {
            let key = (a.min(b), a.max(b));
            if let Some(&(_, m)) = mids.iter().find(|(k, _)| *k == key) {
                return m;
            }
            let (pa, pb) = (v[a as usize], v[b as usize]);
            v.push(unit([0, 1, 2].map(|i| 0.5 * (pa[i] + pb[i]))));
            let m = (v.len() - 1) as u32;
            mids.push((key, m));
            m
        };
        let mut next = Vec::new();
        for [a, b, c] in f {
            let ab = mid(a, b, &mut v);
            let bc = mid(b, c, &mut v);
            let ca = mid(c, a, &mut v);
            next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        f = next;
    }
    let count = f.len();
    SolidMesh {
        id: SolidId(0),
        name: None,
        vertices: v
            .into_iter()
            .map(|p| [0, 1, 2].map(|i| centre[i] + r * p[i]))
            .collect(),
        triangles: f,
        triangle_face: vec![FaceId(0); count],
        exact_volume: None,
    }
}

/// Bounding box of a mesh.
pub(crate) fn bounds(mesh: &SolidMesh) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in &mesh.vertices {
        for a in 0..3 {
            lo[a] = lo[a].min(v[a]);
            hi[a] = hi[a].max(v[a]);
        }
    }
    (lo, hi)
}

/// Total triangle area per face id.
pub(crate) fn face_areas(mesh: &SolidMesh, faces: usize) -> Vec<f64> {
    let mut areas = vec![0.0; faces];
    for (t, f) in mesh.triangles.iter().zip(&mesh.triangle_face) {
        let [a, b, c] = t.map(|i| mesh.vertices[i as usize]);
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let x = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        areas[f.0 as usize] += 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    }
    areas
}

/// A prism: the simple polygon `outline` (counter-clockwise in xy,
/// triangulated by `tris`) extruded from z0 to z1. Faces: 0 bottom,
/// 1 top, 2 sides.
pub(crate) fn extrude(outline: &[[f64; 2]], tris: &[[u32; 3]], z0: f64, z1: f64) -> SolidMesh {
    extrude_faces(outline, tris, z0, z1, &vec![2; outline.len()])
}

/// As [`extrude`], with face `edge_faces[i]` on the side through
/// outline edge i (from point i to point i + 1).
pub(crate) fn extrude_faces(
    outline: &[[f64; 2]],
    tris: &[[u32; 3]],
    z0: f64,
    z1: f64,
    edge_faces: &[u32],
) -> SolidMesh {
    let n = outline.len() as u32;
    let mut vertices: Vec<[f64; 3]> = outline.iter().map(|p| [p[0], p[1], z0]).collect();
    vertices.extend(outline.iter().map(|p| [p[0], p[1], z1]));
    let mut triangles = Vec::new();
    let mut triangle_face = Vec::new();
    for t in tris {
        triangles.push([t[0], t[2], t[1]]);
        triangle_face.push(FaceId(0));
        triangles.push([t[0] + n, t[1] + n, t[2] + n]);
        triangle_face.push(FaceId(1));
    }
    for i in 0..n {
        let j = (i + 1) % n;
        triangles.push([i, j, j + n]);
        triangles.push([i, j + n, i + n]);
        triangle_face.extend([FaceId(edge_faces[i as usize]); 2]);
    }
    SolidMesh {
        id: SolidId(0),
        name: None,
        vertices,
        triangles,
        triangle_face,
        exact_volume: None,
    }
}

/// Fan triangulation of a polygon from vertex `apex` (valid when every
/// vertex is visible from it).
pub(crate) fn fan(n: usize, apex: usize) -> Vec<[u32; 3]> {
    (1..n - 1)
        .map(|k| [apex, (apex + k) % n, (apex + k + 1) % n].map(|v| v as u32))
        .collect()
}

/// Quarter of a plate with a central hole: x in [0, w], y in [0, h],
/// hole radius `r` at the origin (`arc` segments), thickness `t`
/// (z in [0, t]). The outer edges x = w and y = h are split into
/// `outer` segments each. Faces: 0 bottom, 1 top, 2 the y = 0 edge,
/// 3 the x = 0 edge, 4 the hole, then the x = w segments (bottom to
/// top) and the y = h segments (right to left). Returns the mesh, the
/// face count, and each outer segment's (face, midpoint, outward
/// normal).
#[allow(clippy::type_complexity)]
pub(crate) fn quarter_plate_with_hole(
    w: f64,
    h: f64,
    r: f64,
    t: f64,
    arc: usize,
    outer: usize,
) -> (SolidMesh, usize, Vec<(u32, [f64; 2], [f64; 2])>) {
    let mut outline: Vec<[f64; 2]> = vec![[r, 0.0]];
    let mut faces: Vec<u32> = vec![2];
    let mut segments = Vec::new();
    let mut next_face = 5u32;
    for k in 0..outer {
        let y0 = h * k as f64 / outer as f64;
        let y1 = h * (k + 1) as f64 / outer as f64;
        outline.push([w, y0]);
        faces.push(next_face);
        segments.push((next_face, [w, 0.5 * (y0 + y1)], [1.0, 0.0]));
        next_face += 1;
    }
    for k in 0..outer {
        let x0 = w - w * k as f64 / outer as f64;
        let x1 = w - w * (k + 1) as f64 / outer as f64;
        outline.push([x0, h]);
        faces.push(next_face);
        segments.push((next_face, [0.5 * (x0 + x1), h], [0.0, 1.0]));
        next_face += 1;
    }
    outline.push([0.0, h]);
    faces.push(3);
    for k in 0..arc {
        let a = 0.5 * std::f64::consts::PI * (1.0 - k as f64 / arc as f64);
        outline.push(if k == 0 {
            [0.0, r]
        } else {
            [r * a.pcos(), r * a.psin()]
        });
        faces.push(4);
    }
    // The apex (w, h) sees every vertex: the hole is convex toward it.
    let apex = 1 + outer;
    let tris = fan(outline.len(), apex);
    (
        extrude_faces(&outline, &tris, 0.0, t, &faces),
        next_face as usize,
        segments,
    )
}

/// A U-shaped prism 3 × 2 × 1: two prongs separated by a slot of
/// width `slot` (centred at x = 1.5) running from y = 0.5 to the top.
pub(crate) fn u_prism(slot: f64) -> SolidMesh {
    let (a, b) = (1.5 - slot / 2.0, 1.5 + slot / 2.0);
    let outline = [
        [0.0, 0.0],
        [3.0, 0.0],
        [3.0, 2.0],
        [b, 2.0],
        [b, 0.5],
        [a, 0.5],
        [a, 2.0],
        [0.0, 2.0],
    ];
    let tris = [
        [0, 1, 4],
        [1, 2, 3],
        [1, 3, 4],
        [0, 4, 5],
        [0, 5, 6],
        [0, 6, 7],
    ];
    extrude(&outline, &tris, 0.0, 1.0)
}

/// A tube about the z axis: radii `ri` < `ro`, length `length`,
/// `segments` facets around. Faces: 0 inner, 1 outer, 2 ends. Each
/// quad is wound so its normal points away from the material.
pub(crate) fn tube(ri: f64, ro: f64, length: f64, segments: usize) -> SolidMesh {
    let mut vertices = Vec::new();
    for (r, z) in [(ri, 0.0), (ri, length), (ro, 0.0), (ro, length)] {
        for s in 0..segments {
            let t = 2.0 * std::f64::consts::PI * s as f64 / segments as f64;
            vertices.push([r * t.pcos(), r * t.psin(), z]);
        }
    }
    let n = segments as u32;
    let at = |ring: u32, s: usize| ring * n + (s % segments) as u32;
    let mut triangles = Vec::new();
    let mut triangle_face = Vec::new();
    let mut quad = |q: [u32; 4], outward: [f64; 3], face: u32, vertices: &[[f64; 3]]| {
        let p = q.map(|i| vertices[i as usize]);
        let u = [0, 1, 2].map(|k| p[1][k] - p[0][k]);
        let w = [0, 1, 2].map(|k| p[2][k] - p[0][k]);
        let c = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let q = if c[0] * outward[0] + c[1] * outward[1] + c[2] * outward[2] >= 0.0 {
            q
        } else {
            [q[0], q[3], q[2], q[1]]
        };
        triangles.push([q[0], q[1], q[2]]);
        triangles.push([q[0], q[2], q[3]]);
        triangle_face.extend([FaceId(face); 2]);
    };
    for s in 0..segments {
        let mid = 2.0 * std::f64::consts::PI * (s as f64 + 0.5) / segments as f64;
        let radial = [mid.pcos(), mid.psin(), 0.0];
        let inward = radial.map(|c| -c);
        quad(
            [at(0, s), at(0, s + 1), at(1, s + 1), at(1, s)],
            inward,
            0,
            &vertices,
        );
        quad(
            [at(2, s), at(2, s + 1), at(3, s + 1), at(3, s)],
            radial,
            1,
            &vertices,
        );
        quad(
            [at(0, s), at(2, s), at(2, s + 1), at(0, s + 1)],
            [0.0, 0.0, -1.0],
            2,
            &vertices,
        );
        quad(
            [at(1, s), at(3, s), at(3, s + 1), at(1, s + 1)],
            [0.0, 0.0, 1.0],
            2,
            &vertices,
        );
    }
    SolidMesh {
        id: SolidId(0),
        name: None,
        vertices,
        triangles,
        triangle_face,
        exact_volume: None,
    }
}

/// A one-solid analysis model of `mesh` with `faces` faces (areas from
/// the mesh, freeform surfaces).
pub(crate) fn model_of(mesh: SolidMesh, faces: usize) -> crate::model::AnalysisModel {
    use crate::model::{AnalysisModel, FaceInfo, Provenance, RoleSource, Surface, Units};
    let areas = face_areas(&mesh, faces);
    AnalysisModel {
        units: Units::Millimetre,
        chord_tolerance: 0.005,
        faces: (0..faces)
            .map(|f| FaceInfo {
                id: FaceId(f as u32),
                solid: SolidId(0),
                name: None,
                surface: Surface::Freeform,
                role: None,
                role_source: RoleSource::Inferred,
                area: areas[f],
                exact_area: None,
                folded_area: None,
            })
            .collect(),
        solids: vec![mesh],
        concave_edges: Vec::new(),
        provenance: Provenance::Step {
            file_sha256: String::new(),
        },
    }
}

/// A quarter tube about the z axis, angles 0..π/2: radii `ri` < `ro`,
/// length `length`, `segments` facets around the quarter. Faces: 0
/// inner, 1 outer, 2 ends (z = 0 and z = length), 3 the y = 0 side,
/// 4 the x = 0 side. Quads are wound away from the material.
pub(crate) fn quarter_tube(ri: f64, ro: f64, length: f64, segments: usize) -> SolidMesh {
    let n = segments + 1;
    let mut vertices = Vec::new();
    for (r, z) in [(ri, 0.0), (ri, length), (ro, 0.0), (ro, length)] {
        for s in 0..n {
            let t = 0.5 * std::f64::consts::PI * s as f64 / segments as f64;
            let (c, sn) = if s == segments {
                (0.0, 1.0)
            } else {
                (t.pcos(), t.psin())
            };
            vertices.push([r * c, r * sn, z]);
        }
    }
    let at = |ring: usize, s: usize| (ring * n + s) as u32;
    let mut triangles = Vec::new();
    let mut triangle_face = Vec::new();
    let mut quad = |q: [u32; 4], outward: [f64; 3], face: u32, vertices: &[[f64; 3]]| {
        let p = q.map(|i| vertices[i as usize]);
        let u = [0, 1, 2].map(|k| p[1][k] - p[0][k]);
        let w = [0, 1, 2].map(|k| p[2][k] - p[0][k]);
        let c = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let q = if c[0] * outward[0] + c[1] * outward[1] + c[2] * outward[2] >= 0.0 {
            q
        } else {
            [q[0], q[3], q[2], q[1]]
        };
        triangles.push([q[0], q[1], q[2]]);
        triangles.push([q[0], q[2], q[3]]);
        triangle_face.extend([FaceId(face); 2]);
    };
    for s in 0..segments {
        let mid = 0.5 * std::f64::consts::PI * (s as f64 + 0.5) / segments as f64;
        let radial = [mid.pcos(), mid.psin(), 0.0];
        quad(
            [at(0, s), at(0, s + 1), at(1, s + 1), at(1, s)],
            radial.map(|c| -c),
            0,
            &vertices,
        );
        quad(
            [at(2, s), at(2, s + 1), at(3, s + 1), at(3, s)],
            radial,
            1,
            &vertices,
        );
        quad(
            [at(0, s), at(2, s), at(2, s + 1), at(0, s + 1)],
            [0.0, 0.0, -1.0],
            2,
            &vertices,
        );
        quad(
            [at(1, s), at(3, s), at(3, s + 1), at(1, s + 1)],
            [0.0, 0.0, 1.0],
            2,
            &vertices,
        );
    }
    quad(
        [at(0, 0), at(2, 0), at(3, 0), at(1, 0)],
        [0.0, -1.0, 0.0],
        3,
        &vertices,
    );
    let e = segments;
    quad(
        [at(0, e), at(2, e), at(3, e), at(1, e)],
        [-1.0, 0.0, 0.0],
        4,
        &vertices,
    );
    SolidMesh {
        id: SolidId(0),
        name: None,
        vertices,
        triangles,
        triangle_face,
        exact_volume: None,
    }
}
