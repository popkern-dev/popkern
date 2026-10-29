//! The sweep guard of bevel growth (PK-1001): a growing bevel may only
//! meet the topology around it, never a far piece of the body.

use super::{Motion, faces_at, moves};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};
use std::collections::BTreeMap;

/// No vertex that is not on a bevel may lie in the slab its face swept
/// over the region it covered (`swept`: the face's corners at every
/// event): the growth only meets the topology around it.
pub(super) fn sweep_guard(
    start: &BRep,
    out: &BRep,
    motion: &Motion,
    swept: &BTreeMap<FaceHandle, Vec<Point3>>,
    gate: f64,
) -> Result<(), String> {
    for (&face, corners) in swept {
        let (Some(begin), Some((o1, n))) = (motion.plane(face, 0.0), motion.plane(face, 1.0))
        else {
            continue;
        };
        let depth = (o1 - begin.0).dot(&n);
        let hull = planar_hull(corners, &n);
        for (v, vertex) in start.vertices.iter() {
            if faces_at(start, v).contains(&face)
                || moves(start, motion, v)
                || moves(out, motion, v)
                || !out.vertices.contains_key(v)
            {
                continue;
            }
            let p = vertex.position;
            let h = (p - begin.0).dot(&n);
            // Start and end planes included: a far vertex the bevel
            // only touches (its rail through a boss's corner) is met
            // too. Vertices that end on the bevel face move and were
            // skipped above.
            let inside_slab = if depth < 0.0 {
                h <= gate && h >= depth - gate
            } else {
                h >= -gate && h <= depth + gate
            };
            if inside_slab && in_hull(&hull, &p, &n, gate) {
                return Err(format!("the bevel sweeps through a far vertex at {p:?}"));
            }
        }
    }
    Ok(())
}

/// Frame `(u, v)` spanning the plane normal to `n`.
fn frame(n: &Vec3) -> (Vec3, Vec3) {
    let helper = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = n.cross(&helper).normalized().unwrap_or(helper);
    (u, n.cross(&u))
}

/// The convex hull of `points` projected along `n`, counter-clockwise
/// in [`frame`] coordinates (monotone chain).
fn planar_hull(points: &[Point3], n: &Vec3) -> Vec<(f64, f64)> {
    let (u, v) = frame(n);
    let o = Point3::new(0.0, 0.0, 0.0);
    let mut pts: Vec<(f64, f64)> = points
        .iter()
        .map(|p| ((*p - o).dot(&u), (*p - o).dot(&v)))
        .collect();
    pts.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    pts.dedup();
    if pts.len() < 3 {
        return pts;
    }
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    let mut hull: Vec<(f64, f64)> = Vec::new();
    for pass in 0..2 {
        let start = hull.len();
        let iter: Box<dyn Iterator<Item = &(f64, f64)>> = if pass == 0 {
            Box::new(pts.iter())
        } else {
            Box::new(pts.iter().rev())
        };
        for &p in iter {
            while hull.len() >= start + 2
                && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0
            {
                hull.pop();
            }
            hull.push(p);
        }
        hull.pop();
    }
    hull
}

/// Whether `p` projects inside `hull` or onto its boundary (within
/// `gate`): a vertex the grown rail passes through is met.
fn in_hull(hull: &[(f64, f64)], p: &Point3, n: &Vec3, gate: f64) -> bool {
    if hull.len() < 3 {
        return false;
    }
    let (u, v) = frame(n);
    let o = Point3::new(0.0, 0.0, 0.0);
    let q = ((*p - o).dot(&u), (*p - o).dot(&v));
    (0..hull.len()).all(|k| {
        let (a, b) = (hull[k], hull[(k + 1) % hull.len()]);
        let len = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
        len > 0.0 && ((b.0 - a.0) * (q.1 - a.1) - (b.1 - a.1) * (q.0 - a.0)) / len >= -gate
    })
}
