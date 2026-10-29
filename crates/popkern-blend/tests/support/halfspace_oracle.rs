//! Independent exact volume oracle for planar chamfers of CONVEX
//! bodies: the result is the body clipped by each bevel's half-space,
//! computed here by polygon clipping with no kernel code involved.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;

type Poly = Vec<Point3>;

/// Volume of `body` (convex, planar faces) after equal-distance `d`
/// chamfers on `edges` (each given by its endpoints).
pub(crate) fn chamfered_volume(body: &BRep, edges: &[(Point3, Point3)], d: f64) -> f64 {
    let mut faces: Vec<Poly> = body
        .faces
        .keys()
        .map(|f| body.face_outer_loop_positions(f))
        .collect();
    let centre = centroid(&faces.concat());
    let normals: Vec<Vec3> = faces.iter().map(|f| outward(f, &centre)).collect();
    // Wind every polygon counter-clockwise seen from outside, so the
    // signed volume sum and the CCW-built clip caps agree.
    for (f, n) in faces.iter_mut().zip(&normals) {
        let raw = (f[1] - f[0]).cross(&(f[2] - f[0]));
        if raw.dot(n) < 0.0 {
            f.reverse();
        }
    }
    let mut cuts = Vec::new();
    for &(a, b) in edges {
        let on = |f: &Poly, n: &Vec3| [a, b].iter().all(|q| (*q - f[0]).dot(n).abs() < 1e-9);
        let ns: Vec<Vec3> = faces
            .iter()
            .zip(&normals)
            .filter(|(f, n)| on(f, n))
            .map(|(_, n)| *n)
            .collect();
        assert_eq!(ns.len(), 2, "edge {a:?}–{b:?} must bound two faces");
        cuts.push(bevel(a, b, ns[0], ns[1], d));
    }
    for (n, c) in cuts {
        faces = clip_solid(&faces, n, c);
    }
    volume(&faces)
}

fn centroid(points: &[Point3]) -> Point3 {
    let k = points.len() as f64;
    let (x, y, z) = points
        .iter()
        .fold((0.0, 0.0, 0.0), |s, q| (s.0 + q.x, s.1 + q.y, s.2 + q.z));
    Point3::new(x / k, y / k, z / k)
}

fn outward(f: &Poly, centre: &Point3) -> Vec3 {
    let n = (f[1] - f[0])
        .cross(&(f[2] - f[0]))
        .normalized()
        .expect("face normal");
    if n.dot(&(f[0] - *centre)) < 0.0 {
        -n
    } else {
        n
    }
}

/// The bevel half-space `n · x <= c` for edge `a–b` between faces with
/// outward normals `n1`, `n2`: through the two rails, each `d` from the
/// edge along its face, perpendicular to the edge.
fn bevel(a: Point3, b: Point3, n1: Vec3, n2: Vec3, d: f64) -> (Vec3, f64) {
    let t = (b - a).normalized().expect("edge direction");
    let inward = |n: Vec3, other: Vec3| {
        let v = t.cross(&n);
        if v.dot(&other) > 0.0 { -v } else { v }
    };
    let (p1, p2) = (a + inward(n1, n2) * d, a + inward(n2, n1) * d);
    let mut n = t.cross(&(p2 - p1)).normalized().expect("bevel normal");
    if n.dot(&(n1 + n2)) < 0.0 {
        n = -n;
    }
    (n, n.dot(&(p1 - Point3::new(0.0, 0.0, 0.0))))
}

fn side(n: &Vec3, c: f64, q: &Point3) -> f64 {
    n.dot(&(*q - Point3::new(0.0, 0.0, 0.0))) - c
}

fn clip_solid(faces: &[Poly], n: Vec3, c: f64) -> Vec<Poly> {
    let mut out = Vec::new();
    let mut cut: Vec<Point3> = Vec::new();
    for f in faces {
        let mut g = Vec::new();
        for i in 0..f.len() {
            let (a, b) = (f[i], f[(i + 1) % f.len()]);
            let (da, db) = (side(&n, c, &a), side(&n, c, &b));
            if da <= 0.0 {
                g.push(a);
            }
            // A kept vertex ON the plane bounds the cap too (a bevel
            // through an earlier bevel's corner).
            if da.abs() < 1e-12 {
                cut.push(a);
            }
            if (da < 0.0 && db > 0.0) || (da > 0.0 && db < 0.0) {
                let q = a + (b - a) * (da / (da - db));
                g.push(q);
                cut.push(q);
            }
        }
        if g.len() >= 3 {
            out.push(g);
        }
    }
    cut.dedup_by(|a, b| a.distance_to(b) < 1e-12);
    if cut.len() >= 3 {
        let m = centroid(&cut);
        let u = (cut[0] - m).normalized().expect("cut axis");
        let v = n.cross(&u);
        let mut cap: Vec<(f64, Point3)> = cut
            .iter()
            .map(|q| ((*q - m).dot(&v).atan2((*q - m).dot(&u)), *q))
            .collect();
        cap.sort_by(|x, y| x.0.total_cmp(&y.0));
        cap.dedup_by(|x, y| x.1.distance_to(&y.1) < 1e-12);
        out.push(cap.into_iter().map(|(_, q)| q).collect());
    }
    out
}

fn volume(faces: &[Poly]) -> f64 {
    let o = Point3::new(0.0, 0.0, 0.0);
    let mut v = 0.0;
    for f in faces {
        for i in 1..f.len() - 1 {
            v += (f[0] - o).dot(&(f[i] - o).cross(&(f[i + 1] - o)));
        }
    }
    v / 6.0
}
