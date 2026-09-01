use super::*;

fn cubic_row() -> [Point3; 4] {
    [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.3, 0.7, -0.2),
        Point3::new(0.9, 0.4, 0.5),
        Point3::new(1.2, -0.1, 0.1),
    ]
}

fn eval_cubic(row: &[Point3; 4], t: f64) -> Point3 {
    let u = 1.0 - t;
    let b = [u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t];
    let mut acc = [0.0_f64; 3];
    for (w, cp) in b.iter().zip(row) {
        acc[0] += w * cp.x;
        acc[1] += w * cp.y;
        acc[2] += w * cp.z;
    }
    Point3::new(acc[0], acc[1], acc[2])
}

#[test]
fn elevation_is_exact_and_reduction_inverts_it() {
    let cubic = cubic_row();
    let sextic = elevate_cubic_row_to_sextic(cubic);
    for k in 0..=20 {
        let t = k as f64 / 20.0;
        assert!(eval_cubic(&cubic, t).distance_to(&eval_sextic_row(&sextic, t)) < 1e-14);
    }
    let back = reduce_sextic_row_to_cubic(&sextic);
    for (a, b) in cubic.iter().zip(&back) {
        assert!(
            a.distance_to(b) < 1e-12,
            "reduction not exact: {a:?} vs {b:?}"
        );
    }
}

#[test]
fn split_halves_reproduce_the_row() {
    let row = elevate_cubic_row_to_sextic(cubic_row());
    let (first, second) = split_sextic_row(&row);
    for k in 0..=10 {
        let s = k as f64 / 10.0;
        assert!(eval_sextic_row(&first, s).distance_to(&eval_sextic_row(&row, 0.5 * s)) < 1e-14);
        assert!(
            eval_sextic_row(&second, s).distance_to(&eval_sextic_row(&row, 0.5 + 0.5 * s)) < 1e-14
        );
    }
}

fn eval_scalar_cubic(d: [f64; 4], t: f64) -> f64 {
    let u = 1.0 - t;
    u * u * u * d[0] + 3.0 * u * u * t * d[1] + 3.0 * u * t * t * d[2] + t * t * t * d[3]
}

fn eval_vec_sextic(v: &[Vec3; 7], t: f64) -> Vec3 {
    let b = bernstein6(t);
    let mut acc = Vec3::zero();
    for (w, c) in b.iter().zip(v) {
        acc = acc + *c * *w;
    }
    acc
}

#[test]
fn product_reduction_is_exact_on_elevated_cubic_content() {
    let d = [0.1, -0.3, 0.6, 0.2];
    let cubic_v = [
        Vec3::new(0.2, 0.1, 0.0),
        Vec3::new(-0.4, 0.3, 0.2),
        Vec3::new(0.1, -0.2, 0.5),
        Vec3::new(0.3, 0.3, -0.1),
    ];
    let elevated = elevate_cubic_row_to_sextic(core::array::from_fn(|i| {
        Point3::new(cubic_v[i].x, cubic_v[i].y, cubic_v[i].z)
    }));
    let v: [Vec3; 7] =
        core::array::from_fn(|i| Vec3::new(elevated[i].x, elevated[i].y, elevated[i].z));
    let reduced = scalar_cubic_times_sextic_reduced(d, v);
    for k in 0..=40 {
        let t = k as f64 / 40.0;
        let exact = eval_vec_sextic(&v, t) * eval_scalar_cubic(d, t);
        let got = eval_vec_sextic(&reduced, t);
        assert!(
            (exact - got).length() < 1e-13,
            "t={t}: {exact:?} vs {got:?}"
        );
    }
}

#[test]
fn product_reduction_preserves_end_jets_on_genuine_sextic_content() {
    let d = [0.0, 0.0, 0.0, 0.25];
    let v: [Vec3; 7] = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.2, 0.9, 0.1),
        Vec3::new(-0.5, 0.3, 0.4),
        Vec3::new(0.7, -0.6, 0.2),
        Vec3::new(0.1, 0.8, -0.3),
        Vec3::new(-0.2, 0.2, 0.9),
        Vec3::new(0.4, 0.1, 0.0),
    ];
    let reduced = scalar_cubic_times_sextic_reduced(d, v);
    // Degree-9 product coefficients, independently.
    const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    const C9: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    let mut q = [Vec3::zero(); 10];
    for (i, slot) in q.iter_mut().enumerate() {
        for j in 0..4 {
            if i < j || i - j > 6 {
                continue;
            }
            *slot = *slot + v[i - j] * (d[j] * C3[j] * C6[i - j] / C9[i]);
        }
    }
    // Exact Bernstein end jets: value, first, second derivative.
    let jets9 = |a: Vec3, b: Vec3, c: Vec3| (a, (b - a) * 9.0, (c - b * 2.0 + a) * 72.0);
    let jets6 = |a: Vec3, b: Vec3, c: Vec3| (a, (b - a) * 6.0, (c - b * 2.0 + a) * 30.0);
    for (exact, got) in [
        (
            jets9(q[0], q[1], q[2]),
            jets6(reduced[0], reduced[1], reduced[2]),
        ),
        (
            jets9(q[9], q[8], q[7]),
            jets6(reduced[6], reduced[5], reduced[4]),
        ),
    ] {
        assert!((exact.0 - got.0).length() < 1e-14, "value");
        assert!((exact.1 - got.1).length() < 1e-12, "first derivative");
        assert!((exact.2 - got.2).length() < 1e-11, "second derivative");
    }
    // The interior is a genuine least-squares compromise: bounded by
    // the content it cannot represent, never wild.
    for k in 0..=20 {
        let t = k as f64 / 20.0;
        let exact = eval_vec_sextic(&v, t) * eval_scalar_cubic(d, t);
        assert!(
            (exact - eval_vec_sextic(&reduced, t)).length() < 0.05,
            "t={t}"
        );
    }
}
