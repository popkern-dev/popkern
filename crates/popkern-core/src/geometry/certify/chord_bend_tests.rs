use super::super::{C0Breaks, bent_triangle_chord_bound, triangle_chord_bound};
use super::Uv;
use crate::geometry::certify::surface_jet::SurfaceJetSource;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use std::f64::consts::FRAC_PI_2;

fn torus() -> SurfaceGeometry {
    SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::Z,
        major_radius: 3.0,
        minor_radius: 1.0,
    }
}

/// A triangle at the top of the tube with a long ring-direction edge:
/// its chord sags 3.4e-2 inside the tangent plane, the tube only 1.1e-2
/// across it.
fn tube_top() -> [Uv; 3] {
    let v = FRAC_PI_2;
    [(0.0, v - 0.15), (0.3, v - 0.15), (0.15, v + 0.15)]
}

/// Each edge's midpoint displacement taking it to the surface point
/// nearest the chord's midpoint (Gauss–Newton; any bend is valid).
fn feet(s: &SurfaceGeometry, w: [Uv; 3], p: [Point3; 3]) -> [Uv; 3] {
    let mut out = [(0.0, 0.0); 3];
    for (k, d) in out.iter_mut().enumerate() {
        let (a, b) = (k, (k + 1) % 3);
        let m = (0.5 * (w[a].0 + w[b].0), 0.5 * (w[a].1 + w[b].1));
        let c = Point3::new(
            0.5 * (p[a].x + p[b].x),
            0.5 * (p[a].y + p[b].y),
            0.5 * (p[a].z + p[b].z),
        );
        let mut q = m;
        for _ in 0..10 {
            let r = s.point_at(q.0, q.1) - c;
            let (su, sv) = (s.partial_u(q.0, q.1), s.partial_v(q.0, q.1));
            let (g, h, i) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
            let det = g * i - h * h;
            q.0 -= (i * su.dot(&r) - h * sv.dot(&r)) / det;
            q.1 -= (g * sv.dot(&r) - h * su.dot(&r)) / det;
        }
        *d = (q.0 - m.0, q.1 - m.1);
    }
    out
}

/// Edge control points bowing each edge's midpoint by `d`.
fn controls(w: [Uv; 3], d: [Uv; 3]) -> [Option<Uv>; 3] {
    let mut out = [None; 3];
    for (k, c) in out.iter_mut().enumerate() {
        let (a, b) = (w[k], w[(k + 1) % 3]);
        *c = Some((
            0.5 * (a.0 + b.0) + 2.0 * d[k].0,
            0.5 * (a.1 + b.1) + 2.0 * d[k].1,
        ));
    }
    out
}

/// The largest sampled `|L(y) − S(φ(y))|` over the triangle.
fn sampled_bent(s: &SurfaceGeometry, w: [Uv; 3], p: [Point3; 3], d: [Uv; 3]) -> f64 {
    let n = 80;
    let mut worst = 0.0_f64;
    for i in 0..=n {
        for j in 0..=n - i {
            let l = [
                1.0 - (i + j) as f64 / n as f64,
                i as f64 / n as f64,
                j as f64 / n as f64,
            ];
            let mut x = (
                l[0] * w[0].0 + l[1] * w[1].0 + l[2] * w[2].0,
                l[0] * w[0].1 + l[1] * w[1].1 + l[2] * w[2].1,
            );
            for k in 0..3 {
                let bubble = 4.0 * l[k] * l[(k + 1) % 3];
                x.0 += bubble * d[k].0;
                x.1 += bubble * d[k].1;
            }
            let y = Point3::new(
                l[0] * p[0].x + l[1] * p[1].x + l[2] * p[2].x,
                l[0] * p[0].y + l[1] * p[1].y + l[2] * p[2].y,
                l[0] * p[0].z + l[1] * p[1].z + l[2] * p[2].z,
            );
            worst = worst.max(y.distance_to(&s.point_at(x.0, x.1)));
        }
    }
    worst
}

#[test]
fn an_unbent_triangle_is_no_looser_than_the_affine_bound() {
    let s = torus();
    let w = tube_top();
    let p = w.map(|(u, v)| s.point_at(u, v));
    let mut source = SurfaceJetSource::new(&s).expect("source");
    let none = C0Breaks::default();
    let affine = triangle_chord_bound(&mut source, &none, w, p, 1e-3).expect("affine");
    let bent = bent_triangle_chord_bound(&mut source, &none, w, p, [None; 3], 1e-3)
        .expect("bent")
        .expect("applies");
    // The same pairing, bisected one level deeper (BENT_DEPTH): no
    // looser than the affine bound, and still above the sampled
    // deviation.
    let truth = sampled_bent(&s, w, p, [(0.0, 0.0); 3]);
    assert!(
        bent <= affine * (1.0 + 1e-12) && bent >= truth,
        "{bent} vs {affine}, sampled {truth}"
    );
}

#[test]
fn a_tube_top_triangle_bent_onto_its_chords_is_bounded_by_the_tube_sag() {
    let s = torus();
    let w = tube_top();
    let p = w.map(|(u, v)| s.point_at(u, v));
    let mut source = SurfaceJetSource::new(&s).expect("source");
    let none = C0Breaks::default();
    let d = feet(&s, w, p);
    let affine = triangle_chord_bound(&mut source, &none, w, p, 1e-3).expect("affine");
    let bent = bent_triangle_chord_bound(&mut source, &none, w, p, controls(w, d), 1e-3)
        .expect("bent")
        .expect("applies");
    let truth = sampled_bent(&s, w, p, d);
    assert!(affine > 3e-2, "the fixture's in-plane sag: {affine}");
    assert!(bent >= truth, "{bent} under the sampled {truth}");
    assert!(bent < 0.5 * affine && bent < 1.6e-2, "{bent} vs {affine}");
}

#[test]
fn any_bend_bounds_its_own_pairing() {
    // Soundness does not rest on the bend's choice: wild bends (a third
    // of an edge, toward and away from the chord) still bound what
    // they pair.
    let s = torus();
    let w = tube_top();
    let p = w.map(|(u, v)| s.point_at(u, v));
    let mut source = SurfaceJetSource::new(&s).expect("source");
    let none = C0Breaks::default();
    for d in [
        [(0.0, 0.1), (-0.05, 0.02), (0.03, -0.04)],
        [(0.0, -0.1), (0.05, 0.0), (0.0, 0.05)],
        [(0.1, 0.0), (0.0, 0.0), (0.0, 0.0)],
    ] {
        let bent = bent_triangle_chord_bound(&mut source, &none, w, p, controls(w, d), 0.0)
            .expect("bent")
            .expect("applies");
        let truth = sampled_bent(&s, w, p, d);
        assert!(bent >= truth, "{d:?}: {bent} under {truth}");
    }
}

#[test]
fn a_bend_across_a_crease_does_not_apply() {
    let s = torus();
    let w = tube_top();
    let p = w.map(|(u, v)| s.point_at(u, v));
    let mut source = SurfaceJetSource::new(&s).expect("source");
    // A C⁰ line just outside the triangle that the bend reaches over.
    let crease = C0Breaks {
        v: vec![FRAC_PI_2 - 0.16],
        ..C0Breaks::default()
    };
    let d = [(0.0, -0.05), (0.0, 0.0), (0.0, 0.0)];
    let bent =
        bent_triangle_chord_bound(&mut source, &crease, w, p, controls(w, d), 1e-3).expect("bent");
    assert_eq!(bent, None);
}
