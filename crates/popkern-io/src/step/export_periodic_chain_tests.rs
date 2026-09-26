//! Pins for the periodic loop-chain search (PK-558).

use std::f64::consts::{FRAC_PI_2, TAU};

use super::*;

fn meridian(edge: popkern_core::topology::EdgeHandle, u: f64, up: bool) -> ChainUse {
    let (a, b) = if up {
        ((u, -FRAC_PI_2), (u, FRAC_PI_2))
    } else {
        ((u, FRAC_PI_2), (u, -FRAC_PI_2))
    };
    ChainUse {
        edge,
        end_a: a,
        end_b: b,
        a_is_start: Some(true),
        samples: (0..=8)
            .map(|k| {
                let s = k as f64 / 8.0;
                (a.0 + s * (b.0 - a.0), a.1 + s * (b.1 - a.1))
            })
            .collect(),
    }
}

/// A sphere's u is free at its poles (v = ±π/2).
fn sphere_poles(uv: (f64, f64)) -> (bool, bool) {
    ((uv.1.abs() - FRAC_PI_2).abs() < 1e-12, false)
}

fn seam_edge() -> popkern_core::topology::EdgeHandle {
    popkern_core::topology::EdgeHandle::default()
}

#[test]
fn a_seam_loop_through_the_poles_puts_its_two_uses_a_period_apart() {
    // The sphere loop [seam, seam⁻¹], both stored on u = 2π: chaining
    // by cursor alone leaves both uses on one branch (a zero-area face).
    let seam = seam_edge();
    let uses = [meridian(seam, TAU, true), meridian(seam, TAU, false)];
    for want_ccw in [true, false] {
        let (gap, violations, _, winding_ok, assignment) =
            best_assignment(&uses, (Some(TAU), None), want_ccw, &sphere_poles).expect("assignment");
        assert_eq!(gap, 0.0);
        assert_eq!(violations, 0, "seam uses share a branch");
        assert!(winding_ok, "want_ccw={want_ccw}");
        let du = assignment[1].0.0 - assignment[0].0.0;
        assert!((du.abs() - TAU).abs() < 1e-12, "uses {du} apart");
        // The winding picks the side: up at 2π then down at 4π is
        // clockwise in (u, v); counter-clockwise comes down at 0.
        assert_eq!(du < 0.0, want_ccw, "want_ccw={want_ccw}");
    }
}

#[test]
fn without_a_pole_the_chain_keeps_cursor_matching() {
    // An open two-use loop whose junctions are ordinary points: no
    // branch freedom, so the historical closure-first choice stands.
    let seam = seam_edge();
    let uses = [meridian(seam, 1.0, true), meridian(seam, 1.0, false)];
    let (gap, _, _, _, assignment) =
        best_assignment(&uses, (Some(TAU), None), true, &|_| (false, false)).expect("assignment");
    assert_eq!(gap, 0.0);
    assert_eq!(assignment[1].0, (0.0, 0.0));
}
