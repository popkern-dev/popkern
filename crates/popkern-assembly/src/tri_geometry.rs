//! Triangle proximity primitives for the interference check: closest
//! points point–triangle and segment–segment (Ericson, *Real-Time
//! Collision Detection* §5.1), triangle–triangle distance, the
//! segment–triangle crossing test (Möller–Trumbore on a bounded
//! segment), and ray parity against a triangle soup.

use popkern_core::geometry::{Point3, Vec3};

/// A triangle by its corners.
pub type Tri = [Point3; 3];

/// The point of `tri` closest to `p`.
pub fn closest_on_triangle(p: &Point3, tri: &Tri) -> Point3 {
    let [a, b, c] = *tri;
    let ab = b - a;
    let ac = c - a;
    let ap = *p - a;
    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = *p - b;
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return a + ab * v;
    }
    let cp = *p - c;
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return a + ac * w;
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    a + ab * v + ac * w
}

/// Closest points of segments `p1q1` and `p2q2`.
pub fn closest_segments(p1: Point3, q1: Point3, p2: Point3, q2: Point3) -> (Point3, Point3) {
    let d1 = q1 - p1;
    let d2 = q2 - p2;
    let r = p1 - p2;
    let a = d1.dot(&d1);
    let e = d2.dot(&d2);
    let f = d2.dot(&r);
    let eps = 1e-300;
    let (s, t);
    if a <= eps && e <= eps {
        return (p1, p2);
    }
    if a <= eps {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = d1.dot(&r);
        if e <= eps {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = d1.dot(&d2);
            let denom = a * e - b * b;
            let mut s0 = if denom > eps {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let mut t0 = (b * s0 + f) / e;
            if t0 < 0.0 {
                t0 = 0.0;
                s0 = (-c / a).clamp(0.0, 1.0);
            } else if t0 > 1.0 {
                t0 = 1.0;
                s0 = ((b - c) / a).clamp(0.0, 1.0);
            }
            s = s0;
            t = t0;
        }
    }
    (p1 + d1 * s, p2 + d2 * t)
}

/// Whether segment `p..q` crosses triangle `tri` (a transversal hit
/// strictly inside the segment; coplanar contact is not a crossing).
pub fn segment_crosses(p: &Point3, q: &Point3, tri: &Tri) -> bool {
    let dir = *q - *p;
    let e1 = tri[1] - tri[0];
    let e2 = tri[2] - tri[0];
    let h = dir.cross(&e2);
    let det = e1.dot(&h);
    let scale = e1.length() * e2.length() * dir.length();
    if det.abs() <= 1e-12 * scale {
        return false;
    }
    let inv = 1.0 / det;
    let s = *p - tri[0];
    let u = inv * s.dot(&h);
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let qv = s.cross(&e1);
    let v = inv * dir.dot(&qv);
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    let t = inv * e2.dot(&qv);
    t > 1e-12 && t < 1.0 - 1e-12
}

/// Whether two triangles cross (an edge of one passes through the other).
pub fn triangles_cross(a: &Tri, b: &Tri) -> bool {
    (0..3).any(|i| segment_crosses(&a[i], &a[(i + 1) % 3], b))
        || (0..3).any(|i| segment_crosses(&b[i], &b[(i + 1) % 3], a))
}

/// Distance between two triangles and the closest pair (`a` side, `b`
/// side). Crossing triangles are at distance 0.
pub fn triangle_distance(a: &Tri, b: &Tri) -> (f64, Point3, Point3) {
    let mut best = (f64::INFINITY, a[0], b[0]);
    let mut keep = |pa: Point3, pb: Point3| {
        let d = (pa - pb).length();
        if d < best.0 {
            best = (d, pa, pb);
        }
    };
    for p in a {
        keep(*p, closest_on_triangle(p, b));
    }
    for p in b {
        keep(closest_on_triangle(p, a), *p);
    }
    for i in 0..3 {
        for j in 0..3 {
            let (pa, pb) = closest_segments(a[i], a[(i + 1) % 3], b[j], b[(j + 1) % 3]);
            keep(pa, pb);
        }
    }
    if best.0 > 0.0 && triangles_cross(a, b) {
        return (0.0, best.1, best.1);
    }
    best
}

/// Ray `origin + t·dir` (t > 0) against `tri`: the hit parameter.
pub fn ray_hit(origin: &Point3, dir: &Vec3, tri: &Tri) -> Option<f64> {
    let e1 = tri[1] - tri[0];
    let e2 = tri[2] - tri[0];
    let h = dir.cross(&e2);
    let det = e1.dot(&h);
    if det.abs() <= 1e-14 * e1.length() * e2.length() {
        return None;
    }
    let inv = 1.0 / det;
    let s = *origin - tri[0];
    let u = inv * s.dot(&h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(&e1);
    let v = inv * dir.dot(&q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = inv * e2.dot(&q);
    (t > 1e-12).then_some(t)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn distances_between_parallel_and_crossing_triangles() {
        let a = [p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)];
        let b = [p(0.0, 0.0, 0.5), p(1.0, 0.0, 0.5), p(0.0, 1.0, 0.5)];
        assert!((triangle_distance(&a, &b).0 - 0.5).abs() < 1e-15);
        let crossing = [p(0.2, 0.2, -1.0), p(0.2, 0.2, 1.0), p(0.3, 0.4, 0.0)];
        assert_eq!(triangle_distance(&a, &crossing).0, 0.0);
        // Edge to edge, skew.
        let c = [p(2.0, -1.0, 1.0), p(2.0, 1.0, 1.0), p(3.0, 0.0, 1.0)];
        let d = triangle_distance(&a, &c).0;
        assert!((d - 2.0f64.sqrt()).abs() < 1e-12, "{d}");
    }

    #[test]
    fn ray_parity_counts_a_box_from_inside() {
        let tri = [p(-1.0, -1.0, 1.0), p(1.0, -1.0, 1.0), p(0.0, 1.0, 1.0)];
        let up = Vec3::new(0.0, 0.0, 1.0);
        assert!(ray_hit(&p(0.0, 0.0, 0.0), &up, &tri).is_some());
        assert!(ray_hit(&p(0.0, 0.0, 2.0), &up, &tri).is_none());
    }
}
