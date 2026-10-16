//! PK-841 — planar chamfers meeting four or five at a corner build the
//! exact mitre apex when the bevel planes concur.
//!
//! The equal-distance chamfer of every edge of a regular solid is the
//! half-space intersection of the face planes and the bevel planes. On
//! the octahedron (valence 4) and the icosahedron (valence 5) the
//! bevels at a vertex concur by symmetry, so the body has F + E faces
//! and no corner polygon. Non-concurrent bevels (a stretched
//! octahedron, a rectangular pyramid apex) meet along a mitre tree
//! of ridges instead (PK-848). Before PK-841 both declined with an edge
//! `d/√3` (octahedron) off its bevel plane: the per-edge cluster kept
//! the original vertex, and the mitre stage only rebuilt valence 3.
//!
//! The oracle is independent of the chamfer path: every feasible
//! triple-plane point of the half-space system, hulled by
//! `Part::polyhedron`, meshed for its volume.
use popkern_core::geometry::{Point3, Vec3};
use popkern_part::prelude::*;

fn signed_mesh_volume(part: &Part) -> f64 {
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    v6 / 6.0
}

/// `(unit normal, offset)`: the half-space `n · x ≤ h`.
type HalfSpace = (Vec3, f64);

fn to_vec(p: Point3) -> Vec3 {
    Vec3::new(p.x, p.y, p.z)
}

/// Face and equal-distance bevel half-spaces of a convex polyhedral part.
fn chamfer_half_spaces(part: &Part, d: f64) -> Vec<HalfSpace> {
    let b = part.brep();
    let centroid = {
        let n = b.vertices.len() as f64;
        let s = b
            .vertices
            .iter()
            .fold(Vec3::new(0.0, 0.0, 0.0), |s, (_, v)| s + to_vec(v.position));
        s * (1.0 / n)
    };
    // Face planes from their boundary polygons, oriented outward.
    let mut planes = Vec::new();
    let mut face_plane = std::collections::BTreeMap::new();
    for (fh, _) in b.faces.iter() {
        let ps: Vec<Vec3> = b
            .face_outer_loop_positions(fh)
            .into_iter()
            .map(to_vec)
            .collect();
        let mut n = (ps[1] - ps[0]).cross(&(ps[2] - ps[0]));
        n = n * (1.0 / n.dot(&n).sqrt());
        if n.dot(&(ps[0] - centroid)) < 0.0 {
            n = n * -1.0;
        }
        let h = n.dot(&ps[0]);
        planes.push((n, h));
        face_plane.insert(fh, (n, h, ps));
    }
    for (_, e) in b.edges.iter() {
        let he = b.half_edges.get(e.half_edge).expect("half-edge");
        let tw = b.half_edges.get(he.twin).expect("twin");
        let f1 = b.loops.get(he.loop_).expect("loop").face;
        let f2 = b.loops.get(tw.loop_).expect("loop").face;
        let p0 = to_vec(b.vertices.get(he.origin).expect("v").position);
        let p1 = to_vec(b.vertices.get(tw.origin).expect("v").position);
        let axis = (p1 - p0) * (1.0 / (p1 - p0).dot(&(p1 - p0)).sqrt());
        // The point d into each face, perpendicular to the edge.
        let inset = |f| {
            let (n, _, ps): &(Vec3, f64, Vec<Vec3>) = &face_plane[&f];
            let mut u = n.cross(&axis);
            let mid =
                ps.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, p| s + *p) * (1.0 / ps.len() as f64);
            if u.dot(&(mid - p0)) < 0.0 {
                u = u * -1.0;
            }
            p0 + u * d
        };
        let (q1, q2) = (inset(f1), inset(f2));
        let n_sum = face_plane[&f1].0 + face_plane[&f2].0;
        let m = n_sum * (1.0 / n_sum.dot(&n_sum).sqrt());
        assert!(
            m.dot(&(q1 - q2)).abs() < 1e-12,
            "bevel plane holds both rails"
        );
        planes.push((m, m.dot(&q1)));
    }
    planes
}

/// The half-space intersection, hulled independently of any blend code.
fn half_space_oracle(planes: &[HalfSpace]) -> Part {
    let mut points: Vec<Point3> = Vec::new();
    for i in 0..planes.len() {
        for j in i + 1..planes.len() {
            for k in j + 1..planes.len() {
                let (a, b, c) = (planes[i], planes[j], planes[k]);
                let det = a.0.dot(&b.0.cross(&c.0));
                if det.abs() < 1e-9 {
                    continue;
                }
                let x = (b.0.cross(&c.0) * a.1 + c.0.cross(&a.0) * b.1 + a.0.cross(&b.0) * c.1)
                    * (1.0 / det);
                if planes.iter().all(|(n, h)| n.dot(&x) <= h + 1e-10)
                    && !points
                        .iter()
                        .any(|p| (to_vec(*p) - x).dot(&(to_vec(*p) - x)) < 1e-18)
                {
                    points.push(Point3::new(x.x, x.y, x.z));
                }
            }
        }
    }
    Part::polyhedron(&points, Tolerance::default()).expect("oracle hull")
}

fn assert_exact_chamfer(solid: &Part, d: f64, faces: usize, what: &str) {
    let out = solid
        .chamfer(d, |e: EdgeQuery| e)
        .unwrap_or_else(|e| panic!("{what} d={d}: {e}"));
    let oracle = half_space_oracle(&chamfer_half_spaces(solid, d));
    assert_eq!(
        oracle.face_count(),
        faces,
        "{what} d={d}: oracle face count"
    );
    assert_eq!(out.face_count(), faces, "{what} d={d}: face count");
    let issues = out.brep().validate_all();
    assert!(issues.is_empty(), "{what} d={d}: validate_all: {issues:?}");
    assert_eq!(
        out.open_edge_count().expect("mesh"),
        0,
        "{what} d={d}: watertight"
    );
    let b = out.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "{what} d={d}: genus 0");
    let (got, want) = (signed_mesh_volume(&out), signed_mesh_volume(&oracle));
    assert!(
        (got - want).abs() < 1e-9 * want,
        "{what} d={d}: volume {got:.12}, half-space oracle {want:.12}"
    );
    assert_eq!(
        b.vertices.len(),
        oracle.brep().vertices.len(),
        "{what} d={d}: vertices match the half-space oracle (one apex or one mitre tree per corner)"
    );
}

fn octahedron() -> Part {
    let p = |x, y, z| Point3::new(x, y, z);
    Part::polyhedron(
        &[
            p(1.0, 0.0, 0.0),
            p(-1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, -1.0, 0.0),
            p(0.0, 0.0, 1.0),
            p(0.0, 0.0, -1.0),
        ],
        Tolerance::default(),
    )
    .expect("octahedron")
}

fn icosahedron() -> Part {
    let phi = (1.0 + 5.0f64.sqrt()) / 2.0;
    let mut points = Vec::new();
    for a in [-1.0, 1.0] {
        for b in [-phi, phi] {
            points.push(Point3::new(0.0, a, b));
            points.push(Point3::new(a, b, 0.0));
            points.push(Point3::new(b, 0.0, a));
        }
    }
    Part::polyhedron(&points, Tolerance::default()).expect("icosahedron")
}

/// Octahedron stretched along x: the x-axis vertices keep a four-fold
/// symmetric neighbourhood (concurrent bevels), the y- and z-axis
/// vertices do not.
fn stretched_octahedron(a: f64) -> Part {
    let p = |x, y, z| Point3::new(x, y, z);
    Part::polyhedron(
        &[
            p(a, 0.0, 0.0),
            p(-a, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, -1.0, 0.0),
            p(0.0, 0.0, 1.0),
            p(0.0, 0.0, -1.0),
        ],
        Tolerance::default(),
    )
    .expect("stretched octahedron")
}

/// Square pyramid: a concurrent valence-4 apex over valence-3 base
/// corners whose face angles differ (90° base, oblique sides).
fn square_pyramid() -> Part {
    let p = |x, y, z| Point3::new(x, y, z);
    Part::polyhedron(
        &[
            p(-0.5, -0.5, 0.0),
            p(0.5, -0.5, 0.0),
            p(0.5, 0.5, 0.0),
            p(-0.5, 0.5, 0.0),
            p(0.0, 0.0, 1.0),
        ],
        Tolerance::default(),
    )
    .expect("square pyramid")
}

/// Pentagonal trapezohedron (the d10): planar kites need the apex at
/// `z0 (1 + cos 36°) / (1 − cos 36°)` over zig-zag rings at `±z0`.
fn trapezohedron() -> Part {
    let z0 = 0.1;
    let c36 = (36.0f64).to_radians().cos();
    let h = z0 * (1.0 + c36) / (1.0 - c36);
    let mut points = vec![Point3::new(0.0, 0.0, h), Point3::new(0.0, 0.0, -h)];
    for k in 0..10 {
        let a = (36.0 * k as f64).to_radians();
        let z = if k % 2 == 0 { z0 } else { -z0 };
        points.push(Point3::new(a.cos(), a.sin(), z));
    }
    Part::polyhedron(&points, Tolerance::default()).expect("trapezohedron")
}

#[test]
fn every_octahedron_edge_chamfered_is_the_exact_half_space_intersection() {
    for d in [0.05, 0.15, 0.3] {
        assert_exact_chamfer(&octahedron(), d, 8 + 12, "octahedron");
    }
}

#[test]
fn every_icosahedron_edge_chamfered_is_the_exact_half_space_intersection() {
    for d in [0.05, 0.2] {
        assert_exact_chamfer(&icosahedron(), d, 20 + 30, "icosahedron");
    }
}

/// PK-841: oblique valence-3 corners, where a sliver carries rail
/// scaffold vertices between its anchors and the inset corner.
#[test]
fn square_pyramid_base_and_apex_chamfer_exactly() {
    for d in [0.05, 0.1] {
        assert_exact_chamfer(&square_pyramid(), d, 5 + 8, "square pyramid");
    }
}

#[test]
fn every_trapezohedron_edge_chamfered_is_the_exact_half_space_intersection() {
    let solid = trapezohedron();
    assert_eq!(solid.face_count(), 10, "planar kites");
    for d in [0.02, 0.05] {
        assert_exact_chamfer(&solid, d, 10 + 20, "trapezohedron");
    }
}

/// PK-848: on the stretched octahedron the four bevels at a y- or
/// z-axis vertex do not concur: they meet in two triple points joined
/// by a mitre ridge (a tree, not one apex). The exact result is still
/// the half-space intersection, with no corner polygon face.
#[test]
fn non_concurrent_valence_four_bevels_meet_along_a_mitre_ridge() {
    for d in [0.05, 0.1] {
        assert_exact_chamfer(
            &stretched_octahedron(1.4),
            d,
            8 + 12,
            "stretched octahedron",
        );
    }
}

/// Rectangular pyramid: its apex bevels pair off across the long and
/// short sides, so they meet along a ridge too.
fn rectangular_pyramid() -> Part {
    let p = |x, y, z| Point3::new(x, y, z);
    Part::polyhedron(
        &[
            p(-0.6, -0.35, 0.0),
            p(0.6, -0.35, 0.0),
            p(0.6, 0.35, 0.0),
            p(-0.6, 0.35, 0.0),
            p(0.0, 0.0, 0.8),
        ],
        Tolerance::default(),
    )
    .expect("rectangular pyramid")
}

#[test]
fn rectangular_pyramid_apex_chamfers_along_a_mitre_ridge() {
    for d in [0.04, 0.08] {
        assert_exact_chamfer(&rectangular_pyramid(), d, 5 + 8, "rectangular pyramid");
    }
}
