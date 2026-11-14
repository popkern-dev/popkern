use super::*;

/// A polygon on the plane through `origin` spanned by `e1`, `e2` (unit,
/// perpendicular), its loops given in plane coordinates.
fn polygon(origin: Point3, (e1, e2): (Vec3, Vec3), loops: &[&[(f64, f64)]]) -> Polygon {
    let at = |(x, y): (f64, f64)| origin + e1 * x + e2 * y;
    let mut segments = Vec::new();
    let mut flat = Vec::new();
    for l in loops {
        for k in 0..l.len() {
            let (a, b) = (l[k], l[(k + 1) % l.len()]);
            segments.push((at(a), at(b)));
            flat.push((a, b));
        }
    }
    let bx = box_of(&segments).expect("a loop");
    Polygon {
        origin,
        normal: e1.cross(&e2),
        axes: (e1, e2),
        segments,
        flat,
        bx,
    }
}

/// A small deterministic generator in [0, 1).
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn point(&mut self, s: f64) -> Point3 {
        Point3::new(
            s * (self.next() - 0.5),
            s * (self.next() - 0.5),
            s * (self.next() - 0.5),
        )
    }
}

fn frame(g: &mut Lcg) -> (Vec3, Vec3) {
    let a = Vec3::new(g.next() - 0.5, g.next() - 0.5, g.next() - 0.5);
    let b = Vec3::new(g.next() - 0.5, g.next() - 0.5, g.next() - 0.5);
    let e1 = a.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let e2 = e1
        .cross(&b)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (e1, e2.cross(&e1))
}

#[test]
fn segment_distance_is_the_least_over_both_segments() {
    let mut g = Lcg(7);
    for k in 0..400usize {
        let (p0, q0) = (g.point(4.0), g.point(4.0));
        // Some pairs parallel, some degenerate to a point.
        let p1 = if k.is_multiple_of(7) {
            p0
        } else {
            g.point(4.0)
        };
        let q1 = match k % 5 {
            0 => q0 + (p1 - p0) * 0.7,
            1 => q0,
            _ => g.point(4.0),
        };
        let d = segment_distance(p0, p1, q0, q1);
        let n = 200;
        let mut brute = f64::INFINITY;
        for i in 0..=n {
            let p = p0 + (p1 - p0) * (i as f64 / n as f64);
            // The nearest point of a segment to a point is closed form.
            let w = q1 - q0;
            let t = if w.length_squared() > 0.0 {
                ((p - q0).dot(&w) / w.length_squared()).clamp(0.0, 1.0)
            } else {
                0.0
            };
            brute = brute.min(p.distance_to(&(q0 + w * t)));
        }
        assert!(d <= brute + 1e-12, "pair {k}: {d} over the sampled {brute}");
        assert!(
            d >= brute - (p1 - p0).length() / n as f64 - 1e-12,
            "pair {k}: {d} under the sampled {brute}"
        );
    }
}

#[test]
fn a_segment_gap_is_its_least_distance_to_the_region() {
    // A square of side 2 with a square hole of side 1, turned in space.
    let mut g = Lcg(11);
    let outer: &[(f64, f64)] = &[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
    let hole: &[(f64, f64)] = &[(0.5, 0.5), (0.5, 1.5), (1.5, 1.5), (1.5, 0.5)];
    let inside = |x: f64, y: f64| {
        (0.0..=2.0).contains(&x)
            && (0.0..=2.0).contains(&y)
            && !((0.5..1.5).contains(&x) && (0.5..1.5).contains(&y))
    };
    let n = 80;
    let h = 2.0 / n as f64;
    for k in 0..60 {
        let origin = g.point(2.0);
        let axes = frame(&mut g);
        let poly = polygon(origin, axes, &[outer, hole]);
        let (a, b) = (g.point(5.0), g.point(5.0));
        let gap = poly.segment_gap(a, b);
        let mut brute = f64::INFINITY;
        for i in 0..=n {
            for j in 0..=n {
                let (x, y) = (i as f64 * h, j as f64 * h);
                if !inside(x, y) {
                    continue;
                }
                let p = origin + axes.0 * x + axes.1 * y;
                let w = b - a;
                let t = ((p - a).dot(&w) / w.length_squared()).clamp(0.0, 1.0);
                brute = brute.min(p.distance_to(&(a + w * t)));
            }
        }
        assert!(
            gap <= brute + 1e-12,
            "segment {k}: gap {gap} over the sampled {brute}"
        );
        assert!(
            gap >= brute - h * std::f64::consts::SQRT_2 - 1e-12,
            "segment {k}: gap {gap} under the sampled {brute}"
        );
    }
    // Through the hole: no gap of nought though it crosses the plane.
    let poly = polygon(
        Point3::new(0.0, 0.0, 0.0),
        (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
        &[outer, hole],
    );
    let through = poly.segment_gap(Point3::new(1.0, 1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
    assert!((through - 0.5).abs() < 1e-12, "{through}");
    let across = poly.segment_gap(Point3::new(0.25, 1.0, -1.0), Point3::new(0.25, 1.0, 1.0));
    assert_eq!(across, 0.0);
    let above = poly.segment_gap(Point3::new(0.2, 0.2, 0.3), Point3::new(1.8, 0.2, 0.3));
    assert!((above - 0.3).abs() < 1e-12, "{above}");
}

#[test]
fn polygons_cleared_lie_farther_than_the_band() {
    // Triangles about one another: whatever is cleared stays farther
    // than the band at a lattice of both, and both outcomes occur.
    let mut g = Lcg(23);
    let tol = 0.05;
    let (mut cleared, mut kept) = (0, 0);
    let tri: &[(f64, f64)] = &[(0.0, 0.0), (1.0, 0.0), (0.3, 0.8)];
    let n = 24;
    let lattice = |p: &Polygon| -> Vec<Point3> {
        let mut out = Vec::new();
        for i in 0..=n {
            for j in 0..=n - i {
                let (s, t) = (i as f64 / n as f64, j as f64 / n as f64);
                let (x, y) = (s * 1.0 + t * 0.3, t * 0.8);
                out.push(p.origin + p.axes.0 * x + p.axes.1 * y);
            }
        }
        out
    };
    for _ in 0..300 {
        let a = polygon(g.point(0.9), frame(&mut g), &[tri]);
        let b = polygon(g.point(0.9), frame(&mut g), &[tri]);
        if polygons_apart(&a, &b, tol) {
            cleared += 1;
            let (la, lb) = (lattice(&a), lattice(&b));
            for p in &la {
                for q in &lb {
                    assert!(p.distance_to(q) > tol, "cleared pair within the band");
                }
            }
        } else {
            kept += 1;
        }
    }
    assert!(cleared > 30 && kept > 30, "cleared {cleared}, kept {kept}");
    // Parallel and 0.04 apart, one over the other's interior: kept;
    // 0.06 apart: cleared.
    let flat = (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
    let big: &[(f64, f64)] = &[(-2.0, -2.0), (2.0, -2.0), (0.0, 2.0)];
    let low = polygon(Point3::new(0.0, 0.0, 0.0), flat, &[big]);
    for (h, apart) in [(0.04, false), (0.06, true)] {
        let high = polygon(Point3::new(0.0, 0.0, h), flat, &[tri]);
        assert_eq!(polygons_apart(&low, &high, tol), apart, "at {h}");
    }
}

#[test]
fn only_a_planar_face_of_segments_is_a_polygon() {
    use popkern_core::tolerance::Tolerance;
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    for f in cube.faces.keys() {
        let p = polygon_of(&cube, f).expect("a cube face is a polygon");
        assert_eq!(p.segments.len(), 4);
    }
    let rod = BRep::cylinder(1.0, 2.0, tol).expect("rod");
    // Its caps are bounded by circles, its wall is no plane.
    assert!(rod.faces.keys().all(|f| polygon_of(&rod, f).is_none()));
}
