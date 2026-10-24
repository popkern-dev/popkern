//! Composite validation gate for apex relocation.
//!
//! The apex topology rebuild (full-range wire, vertex relocation,
//! base-edge repair — always on since #207 cleanup) needs **a metric
//! apex relocation actually moves** —
//! `validate_all`'s raw error count doesn't (the V↔T coherence it closes
//! doesn't surface there at the kernel's default linear tolerance,
//! and the remaining `InconsistentOrientation` errors are
//! degenerate-normal issues unrelated to V↔T).
//!
//! This module defines `MetricBreakdown`: a per-fixture decomposition
//! of `validate_all` errors into three buckets:
//!
//! * **Structural** — Euler-Poincaré, half-edge pairing/cycles,
//!   vertex consistency, manifold edges, vertex valence.  Must be
//!   zero on any healthy BRep.  Apex relocation doesn't change these.
//! * **Apex-relocation relevant** — `OpenLoop` (loop-closure gap from
//!   curve endpoints disagreeing with vertex positions) and
//!   `ShellGap` (adjacent half-edges' vertex positions diverging
//!   across an edge).  These are the gaps apex V↔T relocation
//!   is designed to close.
//! * **Other geometric** — `InconsistentOrientation` (degenerate
//!   surface normal at parameter midpoint), `CurveSurfaceDeviation`,
//!   `TwistedSurface`, `Precision`, `SelfIntersection`.  Apex relocation
//!   doesn't address these; they need separate work.
//!
//! Combined with the existing `cap_outer_g0` metric
//! (≤ 5e-3 m), the A.7 flip gate becomes:
//!   * structural == 0,
//!   * apex_rebuild_relevant == 0,
//!   * cap_outer_g0 ≤ 5e-3.
//!
//! `Other geometric` count is tracked but not gated — those go
//! through their own work items.

use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence, try_fillet_corner_full};

/// Wrap a fixture build under the apex-relocation flag set.  Returns the
/// post-fillet BRep and the cap face handles.
fn build_with_apex_rebuild(
    is_pyramid: bool,
    radii: &[f64],
) -> Result<popkern_blend::BlendResult, popkern_blend::BlendError> {
    // The apex-relocation passes (wide wire tolerance, face lift,
    // vertex relocation) are unconditional since #207 cleanup.
    // The sphere-fit override remains opt-in (production
    // divergence on variable-radius corners).
    let _g1 = popkern_blend::dispatch_overrides::SPHERE_FIT_BOUNDARY.set_guard(true);

    let (brep, vertex) = if is_pyramid {
        let pyr = BRep::n_gonal_pyramid(radii.len(), 1.0, 0.8, Tolerance::default()).unwrap();
        let vertex = pick_vertex_with_valence(&pyr, radii.len()).unwrap();
        (pyr, vertex)
    } else {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let vertex = cube.vertices.keys().next().unwrap();
        (cube, vertex)
    };
    let edges = incident_edges(&brep, vertex);
    try_fillet_corner_full(&brep, &edges, radii, true, false)
}

/// Decomposition of `validate_all` errors into apex-relocation concern
/// concern buckets.  See module docs for category definitions.
#[derive(Debug, Default)]
pub(crate) struct MetricBreakdown {
    pub structural: usize,
    pub apex_rebuild_relevant: usize,
    pub other_geometric: usize,
    /// Worst gap reported by an `OpenLoop` or `ShellGap` error,
    /// in metres.  0.0 when both are absent.  Useful for ratchets
    /// that watch for regressions in magnitude even when count
    /// stays flat.
    pub worst_apex_rebuild_gap: f64,
}

/// Bucket a `validate_all` error vector by category.
pub(crate) fn breakdown(errors: &[ValidationError]) -> MetricBreakdown {
    let mut out = MetricBreakdown::default();
    for err in errors {
        match err {
            ValidationError::EulerPoincare { .. }
            | ValidationError::DanglingReference { .. }
            | ValidationError::BrokenTwin { .. }
            | ValidationError::BrokenCycle { .. }
            | ValidationError::LoopJunctionMismatch { .. }
            | ValidationError::InconsistentVertex { .. }
            | ValidationError::NonManifoldEdge { .. }
            | ValidationError::LowValence { .. } => out.structural += 1,
            ValidationError::OpenLoop { gap, .. } => {
                out.apex_rebuild_relevant += 1;
                if *gap > out.worst_apex_rebuild_gap {
                    out.worst_apex_rebuild_gap = *gap;
                }
            }
            ValidationError::ShellGap { gap, .. } => {
                out.apex_rebuild_relevant += 1;
                if *gap > out.worst_apex_rebuild_gap {
                    out.worst_apex_rebuild_gap = *gap;
                }
            }
            ValidationError::InconsistentOrientation { .. }
            | ValidationError::WindingDishonest { .. }
            | ValidationError::InnerLoopWinding { .. }
            | ValidationError::MissingPcurves { .. }
            | ValidationError::CanonicalBackward { .. }
            | ValidationError::CurveWindowOutsideSpan { .. }
            | ValidationError::CurveSurfaceDeviation { .. }
            | ValidationError::TwistedSurface { .. }
            | ValidationError::LoopPeriodWinding { .. }
            | ValidationError::PcurveDeviation { .. }
            | ValidationError::Precision { .. }
            | ValidationError::SelfIntersection { .. } => out.other_geometric += 1,
        }
    }
    out
}

/// A.6 cap-kind dispatch lock-in: every uniform-radius fixture must
/// show `apex_rebuild_relevant == 0` under apex relocation.  After §6c-X.A.7's
/// auto-dispatch (suppress sphere-fit substitution on uniform
/// corners; let the original rolling-ball topology stand),
/// the uniform fixtures stay clean. Variable-radius and mixed
/// pyramid fixtures are tracked by diagnostic probes while their
/// repeated-apex support-continuation gaps remain open.
#[test]
fn uniform_fixtures_apex_rebuild_clean() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
        // PK-239: zero. The PK-236 inventory these ceilings carried
        // (pyramid3_uniform = 9 base-corner end-termination chords
        // worst 1.063e-1 + 2 mitred-apex chords; pyramid4_uniform =
        // 12 substituted-boundary cap edges worst 1.5e-2) is closed
        // by the boundary-ring carve: the pyramid apexes and base
        // corners ride the false-sphere-claim gate arm (uniform
        // cylinder strips on mutually oblique axes) and carve exact
        // rings (measured on-surface deviation ≤ 3.8e-14), so every
        // fixture serves validate_all with zero geometric errors.
        max_other_geometric: usize,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: false,
            max_other_geometric: 0,
        },
        Case {
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
            max_other_geometric: 0,
        },
        Case {
            name: "pyramid4_uniform",
            radii: &[0.05, 0.05, 0.05, 0.05],
            is_pyramid: true,
            max_other_geometric: 0,
        },
    ];
    for case in cases {
        let built = build_with_apex_rebuild(case.is_pyramid, case.radii);
        // PK-545: pyramid3_uniform's base corner, which PK-539 declined
        // (its ring folded across the spokes), serves a K-P ring whose
        // cap seams certify G1 over the whole edge
        // (`apex_retaining_welded_contract` pins the certificates).
        let brep = built
            .unwrap_or_else(|e| panic!("{}: must serve, got {e:?}", case.name))
            .brep;
        let errors = brep.validate_all();
        let m = breakdown(&errors);
        eprintln!(
            "[A.6] {}: structural={} pathA={} other={} worst_pathA_gap={:.3e}",
            case.name,
            m.structural,
            m.apex_rebuild_relevant,
            m.other_geometric,
            m.worst_apex_rebuild_gap,
        );
        assert_eq!(
            m.structural,
            0,
            "{}: structural must be 0; errors:\n{}",
            case.name,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
        assert_eq!(
            m.apex_rebuild_relevant,
            0,
            "{}: apex_rebuild_relevant must be 0 under auto-dispatch; errors:\n{}",
            case.name,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
        assert!(
            m.other_geometric <= case.max_other_geometric,
            "{}: other_geometric regression {} (measured ceiling {}):\n{}",
            case.name,
            m.other_geometric,
            case.max_other_geometric,
            errors
                .iter()
                .map(|e| format!("  - {e:?}"))
                .collect::<Vec<_>>()
                .join("\n"),
        );
    }
}

/// PK-548: the uniform pyramid apex (N = 4, r 0.05) and base corner
/// (N = 3, r 0.10) serve K-P rings that certify G1 on every cap seam
/// over the whole edge, and every corner patch claims G1.
///
/// The pyramid4 apex used to fall back to its guide ring. One apex edge
/// is stored base-to-apex, so its strip's ring station, found by a
/// quantized closest-point search, landed 0.0066 off its three
/// neighbours'. That made the ring boundary concave at a hex corner,
/// where no quad chart is regular. With polished stations the strips
/// end symmetrically, and the chord candidate repairs to G1 (measured
/// 2.3e-4 to 5.3e-4 on its cap seams).
#[test]
fn uniform_pyramid_corners_serve_whole_edge_certified_g1_rings() {
    for (fixture, radii) in [
        ("pyramid4_uniform", &[0.05; 4][..]),
        ("pyramid3_uniform", &[0.10; 3][..]),
    ] {
        let result = build_with_apex_rebuild(true, radii)
            .unwrap_or_else(|e| panic!("{fixture}: must serve, got {e:?}"));
        let caps: Vec<_> = result.corner_patches.iter().map(|p| p.face).collect();
        assert_eq!(caps.len(), radii.len(), "{fixture}: one cap patch per side");
        helpers::assert_cap_seams_certify_g1(&result.brep, &caps, fixture);
        for patch in &result.corner_patches {
            assert_eq!(
                patch.achieved_continuity,
                popkern_blend::ContinuityLevel::G1,
                "{fixture}: the claim is the measured G1"
            );
        }
    }
}
