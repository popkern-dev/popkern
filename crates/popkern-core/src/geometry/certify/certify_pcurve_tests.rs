//! Pcurve jets: kernel agreement, dense-sample containment, and the
//! merged affine-polyline branch (PK-542).

use super::super::pcurve_jet::PCURVE_JET_TERMS;
use super::super::*;
use super::{Rng, aligned};
use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::Pcurve;

fn test_pcurves() -> Vec<Pcurve> {
    let nurbs = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.2, 0.5, 0.0),
            Point3::new(0.5, 0.4, 0.0),
            Point3::new(0.7, 0.9, 0.0),
            Point3::new(0.9, 0.6, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        weights: vec![1.0, 0.7, 1.3, 0.9, 1.1, 1.0],
    };
    vec![
        Pcurve::linear((0.1, 0.2), (0.8, -0.4), 2.0, 5.0),
        Pcurve::elliptical((0.5, 0.5), 0.4, 0.2, 0.3, -1.0, 2.5, 0.0, 1.0),
        Pcurve::sinusoidal(0.2, 1.0, 0.3, -0.5, 0.0, 6.0),
        Pcurve {
            geometry: crate::geometry::pcurve::PcurveGeometry::CylinderArcViaRationalQuad {
                theta_a: 0.3,
                theta_b: 2.1,
                v_const: 0.7,
            },
            t_start: 0.0,
            t_end: 1.0,
        },
        Pcurve::from_nurbs(nurbs),
    ]
}

#[test]
fn pcurve_jets_match_the_kernel_and_contain_dense_samples() {
    let mut rng = Rng(3);
    for pcurve in test_pcurves() {
        let source = PcurveJetSource::new(&pcurve).unwrap_or_else(|e| panic!("{e}"));
        let breaks = source.breaks();
        let (t0, t1) = (pcurve.t_start, pcurve.t_end);
        for _ in 0..40 {
            let t = rng.range(t0 + 0.01, t1 - 0.01);
            if breaks.iter().any(|b| (b - t).abs() < 5e-3) {
                continue;
            }
            let jet = source.jet(Iv::point(t)).unwrap_or_else(|e| panic!("{e}"));
            let p = pcurve.point_at(t);
            let at = |t: f64| pcurve.point_at(t);
            let (h1, h2) = (1e-6, 1e-4);
            let d = (at(t + h1) - at(t - h1)) * (1.0 / (2.0 * h1));
            let dd = ((at(t + h2) - at(t)) + (at(t - h2) - at(t))) * (1.0 / (h2 * h2));
            assert!(
                (jet.d[0].u.mid() - p.x).abs() < 1e-12 && (jet.d[0].v.mid() - p.y).abs() < 1e-12
            );
            assert!(
                (jet.d[1].u.mid() - d.x).abs() < 1e-8 && (jet.d[1].v.mid() - d.y).abs() < 1e-8,
                "{:?} d1 {:?} vs {d:?}",
                pcurve.geometry,
                jet.d[1]
            );
            assert!(
                (jet.d[2].u.mid() - dd.x).abs() < 1e-5 && (jet.d[2].v.mid() - dd.y).abs() < 1e-5,
                "d2 {:?} vs {dd:?}",
                jet.d[2]
            );
            // Each higher derivative against differences of the jet's
            // own previous one.
            for k in 3..PCURVE_JET_TERMS {
                let prev_at = |t: f64| {
                    source
                        .jet(Iv::point(t))
                        .map(|j| (j.d[k - 1].u.mid(), j.d[k - 1].v.mid()))
                };
                if let (Ok(p), Ok(m)) = (prev_at(t + h1), prev_at(t - h1)) {
                    let (fu, fv) = ((p.0 - m.0) / (2.0 * h1), (p.1 - m.1) / (2.0 * h1));
                    let scale = 1.0 + fu.abs().max(fv.abs());
                    assert!(
                        (jet.d[k].u.mid() - fu).abs() < 1e-5 * scale
                            && (jet.d[k].v.mid() - fv).abs() < 1e-5 * scale,
                        "{:?} d{k} {:?} vs ({fu}, {fv})",
                        pcurve.geometry,
                        jet.d[k]
                    );
                }
            }
        }
        for trial in 0..20 {
            let width = [0.2, 1e-3][trial % 2] * (t1 - t0);
            let a = rng.range(t0, t1 - width);
            let (a, b) = aligned(a, a + width, &breaks);
            let boxed = source.jet(Iv::new(a, b)).unwrap_or_else(|e| panic!("{e}"));
            for k in 0..=24 {
                let t = a + (b - a) * k as f64 / 24.0;
                let p = pcurve.point_at(t);
                let slack = 1e-13;
                let value = boxed.d[0];
                assert!(
                    value.u.lo - slack <= p.x && p.x <= value.u.hi + slack,
                    "pcurve box [{a}, {b}] misses t = {t}"
                );
                assert!(value.v.lo - slack <= p.y && p.y <= value.v.hi + slack);
                let inner = source.jet(Iv::point(t)).unwrap_or_else(|e| panic!("{e}"));
                // Derivatives past the continuity order are one-sided
                // exactly at a (C²) knot.
                let on_break = breaks.contains(&t);
                for k in 1..PCURVE_JET_TERMS {
                    let (o, i) = (
                        boxed.d[k],
                        if on_break && k >= 3 {
                            boxed.d[k]
                        } else {
                            inner.d[k]
                        },
                    );
                    assert!(o.u.lo <= i.u.hi && i.u.lo <= o.u.hi);
                    assert!(o.v.lo <= i.v.hi && i.v.lo <= o.v.hi);
                }
            }
        }
    }
}

/// PK-542: a degree-1 pcurve that is one affine map up to roundoff
/// (the kernel's fitted cylinder-rail polylines) is enclosed as a
/// single branch, and that enclosure still contains every span's
/// value and slope; a genuine polyline is not merged.
#[test]
fn an_affine_polyline_pcurve_is_one_branch_and_encloses_every_span() {
    let n = 65;
    let control_points: Vec<Point3> = (0..n)
        .map(|k| {
            let s = k as f64 / (n - 1) as f64;
            Point3::new(0.3 + std::f64::consts::TAU * s, 2.0, 0.0)
        })
        .collect();
    let mut knots = vec![0.0];
    knots.extend((0..n).map(|k| k as f64 / (n - 1) as f64));
    knots.push(1.0);
    let affine = Pcurve::from_nurbs(NurbsCurve {
        degree: 1,
        knots: knots.clone(),
        control_points: control_points.clone(),
        weights: vec![1.0; n],
    });
    let source = PcurveJetSource::new(&affine).unwrap_or_else(|e| panic!("{e}"));
    assert!(source.breaks().is_empty(), "merged into one branch");
    let mut rng = Rng(21);
    for _ in 0..50 {
        let a = rng.range(0.0, 0.97);
        let b = a + rng.range(1e-6, 0.03);
        let jet = source.jet(Iv::new(a, b)).unwrap_or_else(|e| panic!("{e}"));
        for k in 0..=20 {
            let t = a + (b - a) * k as f64 / 20.0;
            let p = affine.point_at(t);
            assert!(jet.d[0].contains(p.x, p.y), "value at {t}");
            let h = 1e-7;
            let slope = (affine.point_at((t + h).min(1.0)).x - affine.point_at((t - h).max(0.0)).x)
                / ((t + h).min(1.0) - (t - h).max(0.0));
            assert!((jet.d[1].u.mid() - slope).abs() < 1e-6, "slope at {t}");
        }
    }
    // A genuine polyline keeps its knots.
    let mut bent = control_points;
    bent[32].y += 1e-3;
    let polyline = Pcurve::from_nurbs(NurbsCurve {
        degree: 1,
        knots,
        control_points: bent,
        weights: vec![1.0; n],
    });
    let source = PcurveJetSource::new(&polyline).unwrap_or_else(|e| panic!("{e}"));
    assert!(!source.breaks().is_empty(), "a bent polyline is not merged");
}
