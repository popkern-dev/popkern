//! Irregular convex cones for apex chamfer tests (PK-997, PK-1001):
//! apex off centre, unequal face slopes, so symmetry cannot hide a
//! wrong corner.

use super::chamfer_fixtures::{p, prism, rect};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// The half-space `(x − a) · n ≤ 0` as a box with one face on the plane.
pub(crate) fn half_space(a: Point3, n: Vec3, size: f64) -> BRep {
    let n = n.normalized().expect("normal");
    let helper = if n.z.abs() < 0.9 {
        Vec3::new(0.0, 0.0, 1.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let u = n.cross(&helper).normalized().expect("u");
    let v = n.cross(&u);
    let q = |su: f64, sv: f64| a + u * (su * size) + v * (sv * size);
    let outline = [q(-1.0, -1.0), q(1.0, -1.0), q(1.0, 1.0), q(-1.0, 1.0)];
    BRep::prism(&outline, -n, 2.0 * size, Tolerance::default()).expect("half-space box")
}

pub(crate) const APEX: (f64, f64, f64) = (0.13, -0.07, 1.1);

/// A convex cone with apex [`APEX`] cut from the slab `0 ≤ z ≤ 3` by
/// half-spaces whose outward normals have the given azimuths and
/// elevations (degrees): unequal slopes, apex off centre.
pub(crate) fn cone(faces: &[(f64, f64)]) -> BRep {
    let a = p(APEX.0, APEX.1, APEX.2);
    let mut body = prism(&rect(-3.0, -3.0, 3.0, 3.0), 0.0, 3.0);
    for &(az, el) in faces {
        let (az, el) = (az.to_radians(), el.to_radians());
        let n = Vec3::new(el.cos() * az.cos(), el.cos() * az.sin(), el.sin());
        body = popkern_bool::intersection(body, half_space(a, n, 6.0), Tolerance::default(), None)
            .expect("cone cut");
    }
    body
}

pub(crate) fn cone4() -> BRep {
    cone(&[(5.0, 40.0), (95.0, 52.0), (170.0, 35.0), (265.0, 47.0)])
}

pub(crate) fn cone5() -> BRep {
    cone(&[
        (0.0, 38.0),
        (80.0, 50.0),
        (150.0, 33.0),
        (220.0, 45.0),
        (290.0, 41.0),
    ])
}

pub(crate) fn cone6() -> BRep {
    cone(&[
        (10.0, 42.0),
        (65.0, 36.0),
        (130.0, 50.0),
        (185.0, 39.0),
        (250.0, 44.0),
        (310.0, 34.0),
    ])
}

/// The apex (the highest vertex) and the far ends of its edges, in
/// angular order about it.
pub(crate) fn apex_edges(body: &BRep) -> (Point3, Vec<Point3>) {
    let apex = body
        .vertices
        .values()
        .map(|v| v.position)
        .max_by(|a, b| a.z.total_cmp(&b.z))
        .expect("apex");
    let mut ends: Vec<Point3> = body
        .edges
        .values()
        .filter_map(|e| {
            let he = &body.half_edges[e.half_edge];
            let o = body.vertices[he.origin].position;
            let d = body.vertices[body.half_edges[he.twin].origin].position;
            if o.distance_to(&apex) < 1e-9 {
                Some(d)
            } else if d.distance_to(&apex) < 1e-9 {
                Some(o)
            } else {
                None
            }
        })
        .collect();
    let angle = |q: &Point3| (q.y - apex.y).atan2(q.x - apex.x);
    ends.sort_by(|a, b| angle(a).total_cmp(&angle(b)));
    (apex, ends)
}

/// Every apex selection that leaves an edge sharp, by shape: each
/// single edge, each adjacent pair, a pair around one sharp edge, and
/// all but one.
pub(crate) fn partial_selections(valence: usize) -> Vec<(String, Vec<usize>)> {
    let mut out = Vec::new();
    for i in 0..valence {
        out.push((format!("edge {i}"), vec![i]));
        out.push((format!("pair {i}"), vec![i, (i + 1) % valence]));
    }
    out.push(("around a sharp edge".into(), vec![0, 2]));
    out.push(("all but one".into(), (0..valence - 1).collect()));
    out
}
