//! Tests for the ring-wide cross-field fit (PK-400).

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use super::*;
use popkern_core::geometry::gregory::{bernstein6, split_sextic_row};

fn eval(row: &[Vec3; 7], t: f64) -> Vec3 {
    let b = bernstein6(t);
    let mut acc = Vec3::zero();
    for (w, c) in b.iter().zip(row) {
        acc = acc + *c * *w;
    }
    acc
}

fn eval_points(row: &[Point3; 7], t: f64) -> Point3 {
    popkern_core::geometry::gregory::eval_sextic_row(row, t)
}

/// A planar ring of three straight rim halves per arc on the plane
/// `z = 0` (normal +z), with in-plane degree-6 jet and curvature
/// fields that are polynomial, so the fit must reproduce them.
struct PlanarRing {
    rows: Vec<[[Point3; 7]; 2]>,
    truth: Vec<[([Vec3; 7], [Vec3; 7]); 2]>,
    samples: Vec<[CrossFieldSamples; 2]>,
}

fn planar_ring(n: usize) -> PlanarRing {
    let mut rows = Vec::new();
    let mut truth = Vec::new();
    let mut samples = Vec::new();
    for k in 0..n {
        let angle = |f: f64| core::f64::consts::TAU * (k as f64 + f) / n as f64;
        let point = |f: f64| Point3::new(angle(f).cos(), angle(f).sin(), 0.0);
        let row_between = |f0: f64, f1: f64| -> [Point3; 7] {
            let (a, b) = (point(f0), point(f1));
            core::array::from_fn(|i| {
                let s = i as f64 / 6.0;
                Point3::new(a.x + (b.x - a.x) * s, a.y + (b.y - a.y) * s, 0.0)
            })
        };
        let halves = [row_between(0.0, 0.5), row_between(0.5, 1.0)];
        // Smoothly varying in-plane jet and a small curvature field,
        // defined as ONE degree-6 row over the whole arc and split at
        // m by de Casteljau, so the two halves are C² there.
        let whole = row_between(0.0, 1.0);
        let inward = |i: usize| -> Vec3 {
            let p = whole[i];
            Vec3::new(-p.x, -p.y, 0.0) * 0.02
        };
        let d1_arc: [Vec3; 7] = core::array::from_fn(|i| {
            let g = i as f64 / 6.0;
            inward(i) * (1.0 + 0.1 * g * (1.0 - g))
        });
        let d2_arc: [Vec3; 7] = core::array::from_fn(|i| inward(i) * 0.01);
        let (d1_first, d1_second) = halves_of(&d1_arc);
        let (d2_first, d2_second) = halves_of(&d2_arc);
        let arc_truth = [(d1_first, d2_first), (d1_second, d2_second)];
        let ts: Vec<f64> = (0..25).map(|j| j as f64 / 24.0).collect();
        let arc_samples: Vec<CrossFieldSamples> = arc_truth
            .iter()
            .map(|(d1, d2)| CrossFieldSamples {
                jets: ts.iter().map(|&t| (t, eval(d1, t))).collect(),
                curvatures: ts.iter().map(|&t| (t, eval(d2, t))).collect(),
            })
            .collect();
        rows.push(halves);
        truth.push(arc_truth);
        let [s0, s1] = <[CrossFieldSamples; 2]>::try_from(arc_samples)
            .ok()
            .expect("two");
        samples.push([s0, s1]);
    }
    PlanarRing {
        rows,
        truth,
        samples,
    }
}

fn pieces<'a>(ring: &'a PlanarRing) -> Vec<[RingPiece<'a>; 2]> {
    (0..ring.rows.len())
        .map(|k| {
            core::array::from_fn(|h| RingPiece {
                row0: &ring.rows[k][h],
                samples: &ring.samples[k][h],
                tangency: ring.samples[k][h]
                    .jets
                    .iter()
                    .map(|&(t, j)| (t, Vec3::new(0.0, 0.0, 1.0), j.length()))
                    .collect(),
            })
        })
        .collect()
}

#[test]
fn free_corners_reproduce_polynomial_fields_exactly() {
    let ring = planar_ring(3);
    let corners: Vec<RingCorner> = vec![None; 3];
    let fitted = fit_ring_cross_fields(&pieces(&ring), &corners).expect("solve");
    for (arc, truth) in fitted.iter().zip(&ring.truth) {
        for h in 0..2 {
            for i in 0..7 {
                assert!(
                    (arc[h].0[i] - truth[h].0[i]).length() < 1e-9,
                    "jet coefficient {i} of half {h}"
                );
                assert!(
                    (arc[h].1[i] - truth[h].1[i]).length() < 1e-9,
                    "curvature coefficient {i} of half {h}"
                );
            }
        }
    }
}

#[test]
fn spoke_junctions_are_c2_and_shared_blocks_agree() {
    let ring = planar_ring(4);
    // Corner jets from the neighbouring rims (the neighbour's leg and
    // scaled second difference), as the extraction supplies them.
    let corners: Vec<RingCorner> = (0..4)
        .map(|k| {
            let a = &ring.rows[k][1];
            let b = &ring.rows[(k + 1) % 4][0];
            let jet_a = CornerJet {
                d1: b[1] - b[0],
                d2: ((b[2] - b[1]) - (b[1] - b[0])) * 1.25,
            };
            let jet_b = CornerJet {
                d1: a[5] - a[6],
                d2: ((a[4] - a[5]) - (a[5] - a[6])) * 1.25,
            };
            Some((jet_a, jet_b))
        })
        .collect();
    let fitted = fit_ring_cross_fields(&pieces(&ring), &corners).expect("solve");
    for k in 0..4 {
        // C² at m_k between half 0 (end) and half 1 (start).
        for row in [
            (&fitted[k][0].0, &fitted[k][1].0),
            (&fitted[k][0].1, &fitted[k][1].1),
        ] {
            let (a, b) = row;
            assert!((a[6] - b[0]).length() < 1e-12);
            assert!(((a[6] - a[5]) - (b[1] - b[0])).length() < 1e-12);
            assert!(
                (((a[6] - a[5]) - (a[5] - a[4])) - ((b[2] - b[1]) - (b[1] - b[0]))).length()
                    < 1e-12
            );
        }
        // The two pieces at c_k claim the same corner block.
        let next = (k + 1) % 4;
        let (_, a0) = split_sextic_row(&ring.rows[k][1]);
        let (b0, _) = split_sextic_row(&ring.rows[next][0]);
        let (a1, a2) = split_vec(&fitted[k][1].0, true, &fitted[k][1].1);
        let (b1, b2) = split_vec(&fitted[next][0].0, false, &fitted[next][0].1);
        let depth1 = |r0: Point3, d1: Vec3| r0 + d1 * 0.5;
        let depth2 = |r0: Point3, d1: Vec3, d2: Vec3| r0 + d1 + d2 * 0.2;
        let claims = [
            (depth1(a0[5], a1[5]), depth1(b0[1], b1[1])),
            (depth1(a0[4], a1[4]), depth2(b0[1], b1[1], b2[1])),
            (depth2(a0[5], a1[5], a2[5]), depth1(b0[2], b1[2])),
            (depth2(a0[4], a1[4], a2[4]), depth2(b0[2], b1[2], b2[2])),
        ];
        for (p, q) in claims {
            assert!(
                p.distance_to(&q) < 1e-12,
                "corner block slot disagrees at c_{k}"
            );
        }
        // Corner jets are reproduced exactly.
        let (jet_a, jet_b) = corners[k].expect("shared");
        assert!((fitted[k][1].0[6] - jet_a.d1).length() < 1e-14);
        assert!((fitted[next][0].0[0] - jet_b.d1).length() < 1e-14);
        // Tangency holds at every sampled station.
        for (half, row0) in fitted[k].iter().zip(&ring.rows[k]) {
            for j in 0..=24 {
                let t = j as f64 / 24.0;
                let jet = eval(&half.0, t);
                assert!(jet.z.abs() < 1e-12, "jet leaves the plane at t={t}");
                let _ = eval_points(row0, t);
            }
        }
    }
}

/// De Casteljau halves of a degree-6 vector row at `t = ½`.
fn halves_of(row: &[Vec3; 7]) -> ([Vec3; 7], [Vec3; 7]) {
    let mut levels = [[Vec3::zero(); 7]; 7];
    levels[0] = *row;
    for level in 1..7 {
        for i in 0..(7 - level) {
            levels[level][i] = (levels[level - 1][i] + levels[level - 1][i + 1]) * 0.5;
        }
    }
    (
        core::array::from_fn(|i| levels[i][0]),
        core::array::from_fn(|i| levels[6 - i][i]),
    )
}

fn split_vec(d1: &[Vec3; 7], at_end: bool, d2: &[Vec3; 7]) -> ([Vec3; 7], [Vec3; 7]) {
    let pick = |row: &[Vec3; 7]| -> [Vec3; 7] {
        let (first, second) = halves_of(row);
        if at_end { second } else { first }
    };
    (pick(d1), pick(d2))
}

#[test]
fn incompatible_corner_curvature_stays_out_of_the_jets() {
    // Two strips meeting at c whose sampled curvature fields disagree
    // (a G1 rolling-ball corner): the jets must stay in-plane while
    // the curvature rows absorb the disagreement.
    let mut ring = planar_ring(3);
    for k in 0..3 {
        for (t, c) in ring.samples[k][1].curvatures.iter_mut() {
            *c = *c + Vec3::new(0.0, 0.0, 0.004 * t.powi(2));
        }
    }
    let corners: Vec<RingCorner> = (0..3)
        .map(|k| {
            let a = &ring.rows[k][1];
            let b = &ring.rows[(k + 1) % 3][0];
            Some((
                CornerJet {
                    d1: b[1] - b[0],
                    d2: ((b[2] - b[1]) - (b[1] - b[0])) * 1.25,
                },
                CornerJet {
                    d1: a[5] - a[6],
                    d2: ((a[4] - a[5]) - (a[5] - a[6])) * 1.25,
                },
            ))
        })
        .collect();
    let fitted = fit_ring_cross_fields(&pieces(&ring), &corners).expect("solve");
    for arc in &fitted {
        for half in arc {
            for j in 0..=24 {
                let t = j as f64 / 24.0;
                let jet = eval(&half.0, t);
                // A curvature incompatibility of 20 % of the jet
                // magnitude (twenty times the curvature field itself)
                // leaks ~ε·w² = 2e-5 into the jets; the weights park
                // the disagreement in the curvature rows.
                assert!(
                    jet.z.abs() / jet.length() < 5e-5,
                    "jet tilt {:.2e} at t={t}",
                    jet.z.abs() / jet.length()
                );
            }
        }
    }
}
