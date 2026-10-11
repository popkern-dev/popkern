//! Pinch dissolution by grid candidates cuts exactly where the
//! all-pairs scan cuts (PK-794).

use popkern_core::geometry::Point3;
use popkern_core::topology::EdgeHandle;

use super::super::BoundaryRing;
use popkern_math::PortableMath;

use super::{all_pinch_pairs, dissolve_pinches, dissolve_pinches_among, pinch_candidates};

/// A chain of `lobes` unit squares touching corner to corner, walked
/// along the bottoms and back along the tops: every shared corner is
/// a pinch. `per_side` samples per square side; `shift` rotates the
/// ring's start so pinches fall at different indices.
fn bow_tie_chain(lobes: usize, per_side: usize, shift: usize) -> BoundaryRing {
    let mut uv = Vec::new();
    let side = |from: (f64, f64), to: (f64, f64), uv: &mut Vec<(f64, f64)>| {
        for k in 0..per_side {
            let t = k as f64 / per_side as f64;
            uv.push((from.0 + (to.0 - from.0) * t, from.1 + (to.1 - from.1) * t));
        }
    };
    // Square k spans x in [k, k+1]; its corners touch the neighbours'
    // at (k, 0.5) and (k + 1, 0.5), like diamonds on a string.
    for k in 0..lobes {
        let x = k as f64;
        side((x, 0.5), (x + 0.5, 0.0), &mut uv);
        side((x + 0.5, 0.0), (x + 1.0, 0.5), &mut uv);
    }
    for k in (0..lobes).rev() {
        let x = k as f64;
        side((x + 1.0, 0.5), (x + 0.5, 1.0), &mut uv);
        side((x + 0.5, 1.0), (x, 0.5), &mut uv);
    }
    let len = uv.len();
    uv.rotate_left(shift % len);
    let pos = uv.iter().map(|&(u, v)| Point3::new(u, v, 0.0)).collect();
    let sources = (0..uv.len())
        .map(|k| (EdgeHandle::default(), k as f64, k as f64 + 1.0))
        .collect();
    BoundaryRing { uv, pos, sources }
}

fn rings(r: &[BoundaryRing]) -> Vec<Vec<(u64, u64)>> {
    r.iter()
        .map(|ring| {
            ring.uv
                .iter()
                .map(|&(u, v)| (u.to_bits(), v.to_bits()))
                .collect()
        })
        .collect()
}

#[test]
fn grid_candidates_cut_where_the_all_pairs_scan_cuts() {
    for (lobes, per_side, shift) in [(2, 3, 0), (3, 4, 5), (5, 6, 17), (8, 5, 41), (12, 7, 3)] {
        let ring = || bow_tie_chain(lobes, per_side, shift);
        let fast = dissolve_pinches(ring());
        let reference = dissolve_pinches_among(ring(), all_pinch_pairs);
        assert_eq!(fast.len(), lobes, "{lobes} lobes, shift {shift}");
        assert_eq!(
            rings(&fast),
            rings(&reference),
            "{lobes} lobes, shift {shift}"
        );
    }
}

/// Every pair within `eps_uv` is a candidate, wherever the grid cells
/// fall: pairs are placed at random offsets up to 0.99 eps_uv in
/// random directions, so most straddle a cell boundary.
#[test]
fn every_pair_within_eps_is_a_candidate() {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut unit = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    for trial in 0..40 {
        let n = 60;
        let mut uv: Vec<(f64, f64)> = (0..n).map(|_| (unit() * 3.0, unit() * 2.0)).collect();
        let (lo_u, hi_u, lo_v, hi_v) = uv.iter().fold(
            (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            |(a, b, c, d), &(u, v)| (a.min(u), b.max(u), c.min(v), d.max(v)),
        );
        let eps_uv = (hi_u - lo_u).phypot(hi_v - lo_v) * 1e-6;
        for k in 0..10 {
            let (i, j) = (k * 5, k * 5 + 3 + k % 4);
            let angle = unit() * std::f64::consts::TAU;
            let r = (0.5 + 0.49 * unit()) * eps_uv;
            uv[j] = (uv[i].0 + r * angle.pcos(), uv[i].1 + r * angle.psin());
        }
        let ring = BoundaryRing {
            pos: uv.iter().map(|&(u, v)| Point3::new(u, v, 0.0)).collect(),
            sources: (0..n).map(|_| (EdgeHandle::default(), 0.0, 1.0)).collect(),
            uv,
        };
        let candidates: std::collections::HashSet<(usize, usize)> =
            pinch_candidates(&ring, eps_uv).into_iter().collect();
        for (i, j) in all_pinch_pairs(&ring, eps_uv) {
            let (p, q) = (ring.uv[i], ring.uv[j]);
            if (p.0 - q.0).phypot(p.1 - q.1) <= eps_uv {
                assert!(
                    candidates.contains(&(i, j)),
                    "trial {trial}: pair ({i}, {j}) missed"
                );
            }
        }
    }
}
