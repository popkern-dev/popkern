//! Tests for ring-level pairwise crossing detection.

use super::*;
use crate::geometry::Point3;

/// Build a synthetic pyramid cap ring of `n` bilinear patches obeying the
/// ring convention: apex at `(0, 0, 1)`, spoke `k` a straight segment from
/// the apex to `Q_k = (cos θ_k, sin θ_k, 0)`, patch `k` the bilinear Coons
/// sheet over corners `(0,0) = Q_{k-1}`, `(1,0) = apex`, `(1,1) = Q_k`,
/// `(0,1) = Q_{k-1} + Q_k` (outer corner). Bilinear functions are exactly
/// represented by degree-9 Bernstein grids sampled at Greville abscissae
/// `i/9`, so the grid IS the surface.
fn pyramid_ring(n: usize) -> Vec<RectangularGregoryBiNonic> {
    let apex = Point3::new(0.0, 0.0, 1.0);
    let outer: Vec<Point3> = (0..n)
        .map(|k| {
            let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
            Point3::new(theta.cos(), theta.sin(), 0.0)
        })
        .collect();
    (0..n)
        .map(|k| {
            let q_prev = outer[(k + n - 1) % n];
            let q_next = outer[k];
            let far = Point3::new(q_prev.x + q_next.x, q_prev.y + q_next.y, 0.0);
            let corner = |u: f64, v: f64| -> Point3 {
                let w00 = (1.0 - u) * (1.0 - v);
                let w01 = (1.0 - u) * v;
                let w10 = u * (1.0 - v);
                let w11 = u * v;
                Point3::new(
                    w00 * q_prev.x + w01 * far.x + w10 * apex.x + w11 * q_next.x,
                    w00 * q_prev.y + w01 * far.y + w10 * apex.y + w11 * q_next.y,
                    w00 * q_prev.z + w01 * far.z + w10 * apex.z + w11 * q_next.z,
                )
            };
            let mut patch = RectangularGregoryBiNonic::zeros();
            for i in 0..10 {
                for j in 0..10 {
                    patch.grid[i][j] = corner(i as f64 / 9.0, j as f64 / 9.0);
                }
            }
            patch
        })
        .collect()
}

fn translate(patch: &mut RectangularGregoryBiNonic, dx: f64, dy: f64, dz: f64) {
    for row in patch.grid.iter_mut() {
        for p in row.iter_mut() {
            p.x += dx;
            p.y += dy;
            p.z += dz;
        }
    }
}

#[test]
fn ring_convention_glues_in_the_fixture() {
    // Sanity-pin the fixture itself: spoke k of patch k matches patch
    // (k+1)'s v=0 side reversed, and all (1,0) corners coincide at the
    // apex — the convention the detector documents.
    let ring = pyramid_ring(4);
    for k in 0..4 {
        let a = &ring[k];
        let b = &ring[(k + 1) % 4];
        for s in 0..=8 {
            let t = s as f64 / 8.0;
            let on_a = a.point_at(1.0, t);
            let on_b = b.point_at(1.0 - t, 0.0);
            assert!(
                on_a.distance_to(&on_b) < 1e-12,
                "spoke {k} t={t}: {on_a:?} vs {on_b:?}"
            );
        }
        assert!(
            a.point_at(1.0, 0.0)
                .distance_to(&Point3::new(0.0, 0.0, 1.0))
                < 1e-12
        );
    }
}

#[test]
fn clean_pyramid_rings_report_no_crossing() {
    for n in [3usize, 4, 5] {
        let ring = pyramid_ring(n);
        let crossing = binonic_ring_crossing(&ring, 1.0e-6);
        assert!(
            crossing.is_none(),
            "n={n}: clean ring flagged: {}",
            crossing.map(|c| c.diagnostic_label()).unwrap_or_default()
        );
    }
}

#[test]
fn degenerate_ring_sizes_are_trivially_clean() {
    assert!(binonic_ring_crossing(&[], 1.0e-6).is_none());
    let ring = pyramid_ring(4);
    assert!(binonic_ring_crossing(&ring[..1], 1.0e-6).is_none());
}

#[test]
fn translated_neighbor_patch_is_flagged_as_spoke_neighbor_crossing() {
    // Translating a whole patch keeps it per-patch pristine (the class
    // the per-patch fold gates cannot see) while driving its sheet
    // through the neighbouring sector.
    let mut ring = pyramid_ring(4);
    translate(&mut ring[0], 0.0, 0.4, 0.0);
    // The per-patch gate stays silent on every patch.
    for (i, patch) in ring.iter().enumerate() {
        assert!(
            !GregoryPatch::RectangularBiNonic(Box::new(patch.clone()))
                .self_intersects_with_tolerance(1.0e-6),
            "patch {i} unexpectedly self-intersects"
        );
    }
    let crossing = binonic_ring_crossing(&ring, 1.0e-6)
        .expect("translated neighbour must be flagged as a ring crossing");
    assert_eq!((crossing.patch_a, crossing.patch_b), (0, 1));
    assert_eq!(crossing.relation, RingPairRelation::SpokeNeighbors);
    assert!(
        crossing.separation > crossing.exemption_band,
        "{}",
        crossing.diagnostic_label()
    );
    assert!(crossing.point.x.is_finite());
    assert!(crossing.diagnostic_label().starts_with("ring_crossing"));
}

/// Push patch 0's near-spoke control row (`grid[8][*]`, the same DOF
/// class the exact-boundary optimizer moves) toward the neighbouring
/// sector by `depth`. The `u = 1` boundary row stays fixed, so this
/// bulges the sheet within its own sector rather than crossing.
fn bump_near_spoke_row(ring: &mut [RectangularGregoryBiNonic], depth: f64) {
    for j in 0..10 {
        ring[0].grid[8][j].y += depth;
    }
}

#[test]
fn near_spoke_bulge_is_not_a_crossing() {
    // A near-spoke control-row bump (the DOF class the exact-boundary
    // optimizer moves) bulges the patch WITHIN its own sector: the two
    // sheets still intersect only along the shared spoke (verified by
    // direct SSI during development — the reported curve sits at
    // uv_a.u = 1.0 exactly), so the gate must stay silent at both
    // sub-macro and macroscopic bump depths. Serving rings additionally
    // hold genuine off-spoke hover branches above the 3D band
    // (pyramid3_uniform: 2.38e-5 against band 1e-5, at uv_a.u = 0.9883)
    // — the rail-hover arm exists for those; its integration pin is the
    // fillet corpus (apex_relocation_metrics), which flips if the arm
    // regresses.
    for depth in [1.0e-4, 0.05] {
        let mut ring = pyramid_ring(4);
        bump_near_spoke_row(&mut ring, depth);
        let crossing = binonic_ring_crossing(&ring, 1.0e-6);
        assert!(
            crossing.is_none(),
            "depth {depth}: sector-internal bulge wrongly declined: {}",
            crossing.map(|c| c.diagnostic_label()).unwrap_or_default()
        );
    }
}

#[test]
fn translated_opposite_patch_is_flagged_as_apex_only_crossing() {
    let mut ring = pyramid_ring(4);
    // Patch 2 spans the sector between Q_1 = (0,1,0) and Q_2 = (-1,0,0);
    // pushing it diagonally across the ring centre drives its sheet
    // through patch 0's sector (between Q_3 = (0,-1,0) and Q_0 =
    // (1,0,0)). Pair (0, 1) stays clean, so deterministic scan order
    // reports (0, 2) first.
    translate(&mut ring[2], 1.4, -1.4, 0.0);
    let crossing = binonic_ring_crossing(&ring, 1.0e-6)
        .expect("translated opposite patch must be flagged as a ring crossing");
    assert_eq!((crossing.patch_a, crossing.patch_b), (0, 2));
    assert_eq!(crossing.relation, RingPairRelation::ApexOnly);
    assert!(
        crossing.separation > crossing.exemption_band,
        "{}",
        crossing.diagnostic_label()
    );
}

/// PK-547: the cheap witness agrees with the full gate on every fixture
/// here. It stays silent on clean rings and within-sector bulges, and
/// proves both crossing classes at a refined intersection point beyond
/// the exemption depth.
#[test]
fn the_crossing_witness_agrees_with_the_gate() {
    for n in [3usize, 4, 5] {
        assert!(binonic_ring_crossing_witness(&pyramid_ring(n), 1.0e-6).is_none());
    }
    for depth in [1.0e-4, 0.05] {
        let mut ring = pyramid_ring(4);
        bump_near_spoke_row(&mut ring, depth);
        assert!(binonic_ring_crossing_witness(&ring, 1.0e-6).is_none());
    }
    let mut neighbour = pyramid_ring(4);
    translate(&mut neighbour[0], 0.0, 0.4, 0.0);
    let mut opposite = pyramid_ring(4);
    translate(&mut opposite[2], 1.4, -1.4, 0.0);
    for (ring, pair, relation) in [
        (neighbour, (0, 1), RingPairRelation::SpokeNeighbors),
        (opposite, (0, 2), RingPairRelation::ApexOnly),
    ] {
        let witness = binonic_ring_crossing_witness(&ring, 1.0e-6)
            .unwrap_or_else(|| panic!("the {relation:?} crossing has a witness"));
        assert_eq!((witness.patch_a, witness.patch_b), pair);
        assert_eq!(witness.relation, relation);
        assert!(witness.separation > witness.exemption_band);
        // The witness point lies on both patches.
        let (a, b) = (&ring[pair.0], &ring[pair.1]);
        let on_a = a.point_at(witness.uv_a.0, witness.uv_a.1);
        let on_b = b.point_at(witness.uv_b.0, witness.uv_b.1);
        assert!(
            (on_a - on_b).length() <= 1.0e-6,
            "{}",
            witness.diagnostic_label()
        );
        assert!((on_a - witness.point).length() <= 1.0e-12);
    }
}

/// PK-547: the bounded verdict agrees with the full gate on clean rings
/// (spoke-neighbour contact certifies, apex-only contact settles through
/// the solver), never clears a crossing, and reports an uncertifiable
/// spoke-neighbour pair as undecided without the solver.
#[test]
fn the_bounded_verdict_clears_clean_rings_and_never_a_crossing() {
    for n in [3usize, 4, 5] {
        let verdict = binonic_ring_crossing_bounded(&pyramid_ring(n), 1.0e-6);
        assert!(
            matches!(verdict, BoundedRingCrossing::Clean),
            "n={n}: {verdict:?}"
        );
    }
    let mut neighbour = pyramid_ring(4);
    translate(&mut neighbour[0], 0.0, 0.4, 0.0);
    let verdict = binonic_ring_crossing_bounded(&neighbour, 1.0e-6);
    assert!(
        matches!(verdict, BoundedRingCrossing::Undecided { pair: (0, 1) }),
        "{verdict:?}"
    );
    let mut opposite = pyramid_ring(4);
    translate(&mut opposite[2], 1.4, -1.4, 0.0);
    let verdict = binonic_ring_crossing_bounded(&opposite, 1.0e-6);
    assert!(
        !matches!(verdict, BoundedRingCrossing::Clean),
        "{verdict:?}"
    );
}
