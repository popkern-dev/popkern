//! PK-112 closure gates (formerly the PK-111 follow-up pin): the
//! boolean kernel emits an OCC-exact cube∪sphere (BRepCheck valid,
//! zero free bounds, volume 2.14682 vs the grid-integrated oracle
//! 2.14679), and the tessellator must mesh the FULL-SPHERE wall face
//! — a periodic face whose outer loop is the seam-pair spine (two
//! meridian half-edges through both poles) carrying eight triangular
//! imprint holes.
//!
//! Closure is measured GEOMETRICALLY, at the model's own accuracy
//! class: every face-mesh boundary segment's midpoint must lie within
//! 2× the boolean's linear tolerance of some other mesh segment
//! (own-face matches count — adaptive refinement legitimately leaves
//! hanging nodes ON a coarser cell's edge). Index- or weld-based edge
//! pairing is NOT the gate: face meshes sample shared curves at
//! independent densities, so exact pairing can never close across
//! T-junctions even when the geometry is closed to well under the
//! model tolerance.

use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::{TessellationOptions, tessellate};

const MODEL_LINEAR: f64 = 1.0e-2;
const ORACLE_VOLUME: f64 = 2.14679;

fn cube_union_sphere() -> BRep {
    let t = Tolerance::default().with_linear(MODEL_LINEAR);
    let cube = BRep::unit_cube(t).expect("cube");
    let sphere = BRep::sphere(0.8, t)
        .expect("sphere")
        .transform(&Transform3::translation(0.5, 0.5, 0.5))
        .expect("translate");
    popkern_bool::union(cube, sphere, t, None).expect("union (kernel-exact per OCC)")
}

fn fine_options() -> TessellationOptions {
    TessellationOptions {
        chord_height: 0.001,
        min_subdivisions: 8,
        max_subdivisions: 128,
    }
}

#[test]
fn sphere_union_wall_face_meshes_full_chart() {
    // The regression this gates: the seam-spine outer loop collapsed
    // the wall face's parameter range to a sliver (u ∈ [0, 0.05] of a
    // 2π chart) and the mesh carried 0.06 of the sphere's 8.04 area.
    let out = cube_union_sphere();
    let mesh = tessellate(&out, &fine_options()).expect("tessellate");
    let wall = mesh
        .face_meshes
        .iter()
        .find(|fm| out.faces[fm.face].inner_loops.len() == 8)
        .expect("wall face with eight imprint holes");
    let mut area = 0.0;
    for tri in &wall.triangles {
        let a = wall.vertices[tri[0]].position;
        let b = wall.vertices[tri[1]].position;
        let c = wall.vertices[tri[2]].position;
        area += (b - a).cross(&(c - a)).length() / 2.0;
    }
    assert!(
        area > 7.5 && area < 8.1,
        "wall face must mesh the full periodic chart minus the eight holes \
         (sphere area 8.04, holes ~0.25); got {area:.3}"
    );
}

#[test]
fn sphere_union_mesh_volume_matches_oracle() {
    let out = cube_union_sphere();
    let mesh = tessellate(&out, &fine_options()).expect("tessellate");
    let volume: f64 = mesh
        .triangles
        .iter()
        .map(|tri| {
            let a = mesh.vertices[tri[0]].position;
            let b = mesh.vertices[tri[1]].position;
            let c = mesh.vertices[tri[2]].position;
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum::<f64>()
        .abs();
    assert!(
        (volume - ORACLE_VOLUME).abs() < 5.0e-3,
        "mesh volume should match the grid oracle {ORACLE_VOLUME}; got {volume:.5}"
    );
}

#[test]
fn sphere_union_mesh_is_geometrically_closed_at_model_tolerance() {
    let out = cube_union_sphere();
    let mesh = tessellate(&out, &fine_options()).expect("tessellate");

    // Per-face boundary segments: triangle edges used exactly once
    // after position welding within the face.
    let mut per_face: Vec<Vec<(Point3, Point3)>> = Vec::new();
    for fm in &mesh.face_meshes {
        let key = |p: &Point3| {
            (
                (p.x * 1.0e9).round() as i64,
                (p.y * 1.0e9).round() as i64,
                (p.z * 1.0e9).round() as i64,
            )
        };
        let mut counts: std::collections::HashMap<_, (i64, Point3, Point3)> =
            std::collections::HashMap::new();
        for tri in &fm.triangles {
            for k in 0..3 {
                let (pa, pb) = (
                    fm.vertices[tri[k]].position,
                    fm.vertices[tri[(k + 1) % 3]].position,
                );
                let (ka, kb) = (key(&pa), key(&pb));
                if ka == kb {
                    continue;
                }
                let entry = counts
                    .entry((ka.min(kb), ka.max(kb)))
                    .or_insert((0, pa, pb));
                entry.0 += if ka < kb { 1 } else { -1 };
            }
        }
        per_face.push(
            counts
                .into_values()
                .filter(|&(c, _, _)| c != 0)
                .map(|(_, a, b)| (a, b))
                .collect(),
        );
    }

    let seg_dist = |p: Point3, a: Point3, b: Point3| {
        let ab = b - a;
        let len2 = ab.dot(&ab);
        if len2 <= 0.0 {
            return p.distance_to(&a);
        }
        let s = ((p - a).dot(&ab) / len2).clamp(0.0, 1.0);
        p.distance_to(&(a + ab * s))
    };

    let gate = MODEL_LINEAR * 2.0;
    let mut worst = 0.0f64;
    let mut worst_seg = None;
    for (fi, segs) in per_face.iter().enumerate() {
        for &(a, b) in segs {
            let mid = a + (b - a) * 0.5;
            let mut best = f64::INFINITY;
            for (fj, other) in per_face.iter().enumerate() {
                for &(c, d) in other {
                    if fi == fj
                        && (c.distance_to(&a) + d.distance_to(&b) < 1.0e-12
                            || c.distance_to(&b) + d.distance_to(&a) < 1.0e-12)
                    {
                        continue; // the segment itself
                    }
                    best = best.min(seg_dist(mid, c, d));
                }
            }
            if best > worst {
                worst = best;
                worst_seg = Some((fi, a, b));
            }
        }
    }
    assert!(
        worst <= gate,
        "every mesh boundary segment must sit within 2×linear ({gate}) of a \
         neighboring segment; worst gap {worst:.4} at {worst_seg:?}"
    );
}
