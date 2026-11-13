use super::*;

/// The box holds every dense sample of the curve and is no larger than
/// they are (to rounding), for circles, arcs crossing extremes, arcs
/// between them, tilted ellipses and reversed ranges.
#[test]
fn the_box_is_exact() {
    let c = Point3::new(1.0, -2.0, 3.0);
    let cases = [
        (
            Vec3::new(3f64.sqrt(), 0.0, 0.0),
            Vec3::new(0.0, 3f64.sqrt(), 0.0),
            (0.9, 0.9 + std::f64::consts::TAU),
        ),
        (
            Vec3::new(1.0, 0.3, -0.2),
            Vec3::new(-0.1, 0.8, 0.5),
            (0.2, 2.9),
        ),
        (
            Vec3::new(1.0, 0.3, -0.2),
            Vec3::new(-0.1, 0.8, 0.5),
            (2.9, 0.2),
        ),
        (
            Vec3::new(0.0, 0.0, 2.0),
            Vec3::new(1.0, 0.0, 0.0),
            (0.1, 0.4),
        ),
        (
            Vec3::new(2.0, 1.0, 0.0),
            Vec3::new(0.0, 0.5, 1.5),
            (-7.0, 9.0),
        ),
    ];
    for (a, b, range) in cases {
        let bx = trig_curve_box(c, a, b, range);
        let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
        let mut seen = Aabb {
            min: c + a * lo.cos() + b * lo.sin(),
            max: c + a * lo.cos() + b * lo.sin(),
        };
        for i in 0..=100_000 {
            let t = lo + (hi - lo) * f64::from(i) / 100_000.0;
            let p = c + a * t.cos() + b * t.sin();
            for (v, l, h) in [
                (p.x, bx.min.x, bx.max.x),
                (p.y, bx.min.y, bx.max.y),
                (p.z, bx.min.z, bx.max.z),
            ] {
                assert!(
                    v >= l - 1e-12 && v <= h + 1e-12,
                    "{v} outside [{l}, {h}] for {range:?}"
                );
            }
            seen = Aabb {
                min: Point3::new(
                    seen.min.x.min(p.x),
                    seen.min.y.min(p.y),
                    seen.min.z.min(p.z),
                ),
                max: Point3::new(
                    seen.max.x.max(p.x),
                    seen.max.y.max(p.y),
                    seen.max.z.max(p.z),
                ),
            };
        }
        assert!(
            bx.min.distance_to(&seen.min) < 1e-8 && bx.max.distance_to(&seen.max) < 1e-8,
            "{range:?}: {bx:?} against {seen:?}"
        );
    }
}
