//! Soundness of the certify substrate: dense samples always land in
//! the enclosures, and point enclosures agree with the kernel's own
//! evaluators (formula correctness).

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::polyhedral_spline_patch::{PolyhedralSplinePatch, PolyhedralSplineProvenance};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use std::ops::Div;

/// Deterministic xorshift in [0, 1).
pub(super) struct Rng(pub(super) u64);

impl Rng {
    pub(super) fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    pub(super) fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next()
    }
}

#[test]
fn interval_arithmetic_encloses_random_point_operations() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for _ in 0..20_000 {
        let (a0, a1) = (rng.range(-3.0, 3.0), rng.range(-3.0, 3.0));
        let (b0, b1) = (rng.range(-3.0, 3.0), rng.range(-3.0, 3.0));
        let (a, b) = (Iv::new(a0, a1), Iv::new(b0, b1));
        let (x, y) = (rng.range(a.lo, a.hi), rng.range(b.lo, b.hi));
        assert!((a + b).contains(x + y));
        assert!((a - b).contains(x - y));
        assert!((a * b).contains(x * y));
        assert!(a.sqr().contains(x * x));
        assert!(a.abs().contains(x.abs()));
        if !b.contains_zero() {
            assert!(a.div(b).contains(x / y));
        }
        if a.lo >= 0.0 {
            assert!(a.sqrt().contains(x.sqrt()));
        }
        assert!(a.atan().contains(x.atan()));
        let wide = Iv::new(a0 * 4.0, a1 * 4.0);
        let z = rng.range(wide.lo, wide.hi);
        assert!(wide.cos().contains(z.cos()), "cos {wide:?} at {z}");
        assert!(wide.sin().contains(z.sin()), "sin {wide:?} at {z}");
    }
    assert_eq!(Iv::new(1.0, 2.0).div(Iv::new(-1.0, 1.0)), Iv::ENTIRE);
    // Extrema inside an interval are reached.
    assert_eq!(Iv::new(-0.1, 0.1).cos().hi, 1.0);
    assert_eq!(Iv::new(3.0, 3.3).cos().lo, -1.0);
    assert_eq!(Iv::new(1.5, 1.6).sin().hi, 1.0);
}

#[test]
fn taylor_arithmetic_matches_closed_form_series() {
    // Series about s = 0.5 of the identity: 0.5 + δ.
    let mut c = [Iv::point(0.0); 5];
    c[0] = Iv::point(0.5);
    c[1] = Iv::point(1.0);
    let s = Taylor::new(c);
    let one = Taylor::<5>::constant(Iv::point(1.0));
    // 1/(1 + s): coefficients (−1)^k / 1.5^{k+1}.
    let q = one / (one + s);
    for k in 0..5 {
        let exact = (-1.0_f64).powi(k as i32) / 1.5_f64.powi(k as i32 + 1);
        assert!((q.c[k].mid() - exact).abs() < 1e-15, "1/(1+s) c{k}");
    }
    // √(1 + s): binomial series coefficients C(1/2, k) · 1.5^{1/2 − k}.
    let r = (one + s).sqrt();
    let mut binom = 1.0;
    for k in 0..5 {
        let exact = binom * 1.5_f64.powf(0.5 - k as f64);
        assert!((r.c[k].mid() - exact).abs() < 1e-14, "sqrt c{k}");
        binom *= (0.5 - k as f64) / (k as f64 + 1.0);
    }
    // atan(s): derivatives of atan at 0.5 via 1/(1 + s²) integrated.
    let a = s.atan();
    assert!((a.c[0].mid() - 0.5_f64.atan()).abs() < 1e-15);
    assert!((a.c[1].mid() - 1.0 / 1.25).abs() < 1e-15);
    assert!((a.c[2].mid() - (-0.5 / 1.25_f64.powi(2))).abs() < 1e-15);
    // Composition: exp-free check (1 + Δ)² composed.
    let square = (one + s).sqr();
    assert!((square.c[0].mid() - 2.25).abs() < 1e-15 && (square.c[2].mid() - 1.0).abs() < 1e-15);
}

fn rational_nurbs() -> NurbsSurface {
    // Degree (3, 2), 5 × 4 controls, two u-spans and two v-spans.
    let mut rng = Rng(42);
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for i in 0..5 {
        for j in 0..4 {
            control_points.push(Point3::new(
                i as f64 + rng.range(-0.2, 0.2),
                j as f64 + rng.range(-0.2, 0.2),
                rng.range(-1.0, 1.0),
            ));
            weights.push(rng.range(0.6, 1.8));
        }
    }
    NurbsSurface {
        degree_u: 3,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 0.55, 1.0, 1.0, 1.0],
        control_points,
        num_u: 5,
        num_v: 4,
        weights,
    }
}

pub(super) fn binonic() -> SurfaceGeometry {
    let mut rng = Rng(7);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, p) in row.iter_mut().enumerate() {
            *p = Point3::new(i as f64 / 9.0, j as f64 / 9.0, 0.3 * rng.range(-1.0, 1.0));
        }
    }
    SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(Box::new(
        RectangularGregoryBiNonic { grid },
    )))
}

fn polyhedral() -> SurfaceGeometry {
    let mut rng = Rng(11);
    let mut net = [[Point3::new(0.0, 0.0, 0.0); 5]; 5];
    for (j, row) in net.iter_mut().enumerate() {
        for (i, p) in row.iter_mut().enumerate() {
            *p = Point3::new(i as f64, j as f64, rng.range(-1.0, 1.0));
        }
    }
    SurfaceGeometry::PolyhedralSpline(Box::new(PolyhedralSplinePatch::new(
        net,
        PolyhedralSplineProvenance {
            valence: 5,
            sector: 0,
        },
    )))
}

type SurfaceCase = (SurfaceGeometry, (f64, f64), (f64, f64));

fn test_surfaces() -> Vec<SurfaceCase> {
    let axis = Vec3::new(0.3, -0.5, 0.8)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    vec![
        (
            SurfaceGeometry::Plane {
                origin: Point3::new(1.0, 2.0, 3.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            (-2.0, 2.0),
            (-1.0, 3.0),
        ),
        (
            SurfaceGeometry::Cylinder {
                origin: Point3::new(0.5, 0.0, -1.0),
                axis,
                radius: 2.5,
            },
            (-1.0, 7.0),
            (-2.0, 2.0),
        ),
        (
            SurfaceGeometry::Cone {
                apex: Point3::new(0.0, 1.0, 0.0),
                axis,
                half_angle: 0.4,
            },
            (0.0, 6.3),
            (0.2, 3.0),
        ),
        (
            SurfaceGeometry::Sphere {
                center: Point3::new(1.0, 1.0, 1.0),
                radius: 1.7,
                axis: crate::geometry::Vec3::Z,
            },
            (-3.0, 3.0),
            (0.1, 3.0),
        ),
        (
            SurfaceGeometry::Torus {
                center: Point3::new(0.0, 0.0, 0.0),
                axis,
                major_radius: 3.0,
                minor_radius: 0.8,
            },
            (0.0, 6.3),
            (-3.2, 3.2),
        ),
        (
            SurfaceGeometry::Nurbs(rational_nurbs()),
            (0.0, 1.0),
            (0.0, 1.0),
        ),
        (binonic(), (0.0, 1.0), (0.0, 1.0)),
        (polyhedral(), (0.0, 1.0), (0.0, 1.0)),
    ]
}

/// Split `[lo, hi]` so no piece straddles a break.
pub(super) fn aligned(lo: f64, hi: f64, breaks: &[f64]) -> (f64, f64) {
    let mut hi = hi;
    for &b in breaks {
        if b > lo && b < hi {
            hi = b;
        }
    }
    (lo, hi)
}

/// Reference partials up to order 2 from the kernel's evaluators
/// (second order by central differences of its first partials).
fn reference(s: &SurfaceGeometry, u: f64, v: f64) -> Vec<((usize, usize), Vec3, f64)> {
    let h = 1e-6;
    let (su, sv) = (|u, v| s.partial_u(u, v), |u, v| s.partial_v(u, v));
    let d1 = |f: &dyn Fn(f64, f64) -> Vec3, du: f64, dv: f64| {
        (f(u + du, v + dv) - f(u - du, v - dv)) * (1.0 / (2.0 * h))
    };
    let p = s.point_at(u, v);
    vec![
        ((0, 0), Vec3::new(p.x, p.y, p.z), 1e-12),
        ((1, 0), su(u, v), 1e-9),
        ((0, 1), sv(u, v), 1e-9),
        ((2, 0), d1(&su, h, 0.0), 1e-6),
        ((1, 1), d1(&su, 0.0, h), 1e-6),
        ((0, 2), d1(&sv, 0.0, h), 1e-6),
    ]
}

/// Third partials against central differences of the jet's own
/// second partials (the Leibniz recursion's self-consistency).
fn check_third_order(source: &mut SurfaceJetSource, name: &str, u: f64, v: f64) {
    let h = 1e-6;
    let at = |source: &mut SurfaceJetSource, u: f64, v: f64| {
        source
            .jet(Iv::point(u), Iv::point(v))
            .unwrap_or_else(|e| panic!("{e}"))
    };
    let (up, um) = (at(source, u + h, v), at(source, u - h, v));
    let (vp, vm) = (at(source, u, v + h), at(source, u, v - h));
    let jet = at(source, u, v);
    let diff = |p: &SurfaceJet, m: &SurfaceJet, a: usize, b: usize| {
        (p.partial(a, b).mid() - m.partial(a, b).mid()) * (1.0 / (2.0 * h))
    };
    for ((a, b), expected) in [
        ((3, 0), diff(&up, &um, 2, 0)),
        ((2, 1), diff(&vp, &vm, 2, 0)),
        ((1, 2), diff(&up, &um, 0, 2)),
        ((0, 3), diff(&vp, &vm, 0, 2)),
    ] {
        let err = (jet.partial(a, b).mid() - expected).length();
        assert!(
            err <= 1e-5 * (1.0 + expected.length()),
            "{name} ∂({a},{b}) at ({u}, {v}): {err:e}"
        );
    }
}

#[test]
fn surface_point_jets_match_the_kernel_evaluators() {
    let mut rng = Rng(5);
    for (surface, ur, vr) in test_surfaces() {
        let mut source = SurfaceJetSource::new(&surface).unwrap_or_else(|e| panic!("{e}"));
        let (bu, bv) = source.breaks();
        for _ in 0..40 {
            let (u, v) = (
                rng.range(ur.0 + 0.01, ur.1 - 0.01),
                rng.range(vr.0 + 0.01, vr.1 - 0.01),
            );
            // Stay clear of knots so the differences see one span.
            if bu.iter().any(|b| (b - u).abs() < 5e-3) || bv.iter().any(|b| (b - v).abs() < 5e-3) {
                continue;
            }
            let jet = source
                .jet(Iv::point(u), Iv::point(v))
                .unwrap_or_else(|e| panic!("{} at ({u}, {v}): {e}", surface.kind_name()));
            for ((a, b), expected, tol) in reference(&surface, u, v) {
                let got = jet.partial(a, b);
                let scale = 1.0 + expected.length();
                let err = (got.mid() - expected).length();
                assert!(
                    err <= tol * scale,
                    "{} ∂({a},{b}) at ({u}, {v}): {err:e} (enclosure {got:?}, kernel {expected:?})",
                    surface.kind_name()
                );
                assert!(
                    got.mag() - got.mid().length() < 1e-9 * scale,
                    "point jet too wide"
                );
            }
            check_third_order(&mut source, surface.kind_name(), u, v);
        }
    }
}

#[test]
fn surface_box_jets_contain_dense_point_jets() {
    let mut rng = Rng(99);
    for (surface, ur, vr) in test_surfaces() {
        let mut source = SurfaceJetSource::new(&surface).unwrap_or_else(|e| panic!("{e}"));
        let (bu, bv) = source.breaks();
        for trial in 0..30 {
            let width = [0.3, 0.05, 1e-4][trial % 3];
            let u0 = rng.range(ur.0, ur.1 - width);
            let v0 = rng.range(vr.0, vr.1 - width);
            let (u0, u1) = aligned(u0, u0 + width, &bu);
            // Alternate full boxes with iso-parametric (thin) boxes.
            let (v0, v1) = if trial % 2 == 0 {
                aligned(v0, v0 + width, &bv)
            } else {
                (v0, v0)
            };
            let boxed = source
                .jet(Iv::new(u0, u1), Iv::new(v0, v1))
                .unwrap_or_else(|e| panic!("{} box: {e}", surface.kind_name()));
            for k in 0..=12 {
                for l in 0..=12 {
                    let u = u0 + (u1 - u0) * k as f64 / 12.0;
                    let v = v0 + (v1 - v0) * l as f64 / 12.0;
                    let point = source
                        .jet(Iv::point(u), Iv::point(v))
                        .unwrap_or_else(|e| panic!("{e}"));
                    let (on_u_break, on_v_break) = (bu.contains(&u), bv.contains(&v));
                    for &(a, b) in surface_jet::JET_ORDERS.iter() {
                        // Partials of order ≥ 2 across a knot may be
                        // one-sided (the fixture is C¹ in v, C² in u):
                        // the box encloses its own span's side.
                        if (on_u_break && a >= 2) || (on_v_break && b >= 2) {
                            continue;
                        }
                        let (outer, inner) = (boxed.partial(a, b), point.partial(a, b));
                        let inside = |o: Iv, i: Iv| o.lo <= i.hi && i.lo <= o.hi;
                        assert!(
                            inside(outer.x, inner.x)
                                && inside(outer.y, inner.y)
                                && inside(outer.z, inner.z),
                            "{} ∂({a},{b}) box [{u0},{u1}]×[{v0},{v1}] misses ({u},{v})",
                            surface.kind_name()
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn nurbs_boxes_straddling_a_knot_are_declined() {
    let mut source = SurfaceJetSource::new(&SurfaceGeometry::Nurbs(rational_nurbs()))
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(
        source.jet(Iv::new(0.3, 0.5), Iv::point(0.2)),
        Err(CertifyDecline::SpanStraddle)
    );
    // Closed spans: a box ending exactly on the knot is inside.
    assert!(source.jet(Iv::new(0.3, 0.4), Iv::point(0.2)).is_ok());
}

#[path = "certify_pcurve_tests.rs"]
mod pcurve;

#[path = "certify_gregory_tests.rs"]
mod gregory;
