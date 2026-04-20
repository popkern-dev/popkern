//! §7 continuity-downgrade plumbing.
//!
//! Collected in one module so the warning-emission and fallback-
//! tagging logic don't bloat `lib.rs` and `geometry_compute.rs`
//! past the 500-line structural cap. See
//! [`crate::types::ContinuityLevel`] and
//! [`crate::error::BlendWarning::ContinuityDowngraded`] for the
//! public surface.

use popkern_core::topology::EdgeHandle;

use crate::construct::EdgeGeometry;
use crate::error::BlendWarning;
use crate::types::{ContinuityLevel, CornerPatchInfo, FilletOptions};

/// Append one `ContinuityDowngraded` warning per edge whose
/// `EdgeGeometry::continuity_downgrade` is populated.
pub(crate) fn emit_downgrade_warnings(
    edge_geometries: &[EdgeGeometry],
    options: &FilletOptions,
    warnings: &mut Vec<(EdgeHandle, BlendWarning)>,
) {
    for geom in edge_geometries {
        if let Some((achieved, reason)) = &geom.continuity_downgrade {
            let edge = geom.analysis.edge;
            warnings.push((
                edge,
                BlendWarning::ContinuityDowngraded {
                    edge,
                    requested: options.continuity,
                    achieved: *achieved,
                    reason: reason.clone(),
                },
            ));
        }
    }
}

/// Promote a walker/sweep `downgrade_reason` string into the typed
/// `(achieved, reason)` pair stored on `EdgeGeometry`. Returns
/// `None` when no downgrade occurred. The achieved level is always
/// `G1` in the current pipeline — the fallback paths land on arc
/// sweeps / analytic torus output.
pub(crate) fn to_field(downgrade_reason: Option<String>) -> Option<(ContinuityLevel, String)> {
    downgrade_reason.map(|r| (ContinuityLevel::G1, r))
}

/// Append one `CornerContinuityDowngraded` warning **per vertex**
/// (not per patch) when the corner falls short of the caller's
/// requested continuity — either (a) achieved < requested, or
/// (b) `n_sides >= 5` (a quality caveat at high valence).
///
/// **R8 C7 dedup.** When the three-quad dispatch is on, a single
/// vertex can contribute up to 3 `CornerPatchInfo` entries (one
/// per sub-patch). Emitting three warnings about the same
/// geometric corner is duplicative and noisy. We now aggregate:
/// the first entry per vertex drives the warning, subsequent
/// entries for the same vertex are skipped. When the sub-patches
/// report different achieved levels (shouldn't happen today —
/// `apply_three_quad_decomposition_to_blend` labels all three
/// G1 — but might in future with twist extraction), the
/// aggregated report uses the MIN achieved level so the warning
/// is honest about the weakest sub-patch.
///
/// `fallback_edge` is used only as the outer tuple key so the
/// warning rides on the existing
/// `Vec<(EdgeHandle, BlendWarning)>` channel; the authoritative
/// location is the `vertex` field inside the variant body.
///
/// When both conditions hold (N≥5 *and* achieved < requested),
/// a single warning is emitted — the N≥5 rationale subsumes the
/// downgrade rationale at high valence.
pub(crate) fn emit_corner_downgrade_warnings(
    corner_patches: &[CornerPatchInfo],
    options: &FilletOptions,
    fallback_edge: EdgeHandle,
    warnings: &mut Vec<(EdgeHandle, BlendWarning)>,
) {
    use std::collections::BTreeMap;

    let requested = options.continuity;

    // Aggregate per-vertex: track min achieved continuity + the
    // first patch's (n_sides, method) for the warning message.
    // BTreeMap gives deterministic iteration — important for
    // stable warning order across runs.
    struct Agg<'a> {
        first: &'a CornerPatchInfo,
        min_achieved: crate::types::ContinuityLevel,
    }
    let mut by_vertex: BTreeMap<_, Agg<'_>> = BTreeMap::new();
    for info in corner_patches {
        by_vertex
            .entry(info.vertex)
            .and_modify(|agg| {
                if info.achieved_continuity < agg.min_achieved {
                    agg.min_achieved = info.achieved_continuity;
                }
            })
            .or_insert(Agg {
                first: info,
                min_achieved: info.achieved_continuity,
            });
    }

    for agg in by_vertex.into_values() {
        let info = agg.first;
        if info.n_sides >= 5 {
            warnings.push((
                fallback_edge,
                BlendWarning::CornerContinuityDowngraded {
                    vertex: info.vertex,
                    n_sides: info.n_sides,
                    requested,
                    achieved: agg.min_achieved,
                    reason: format!(
                        "N = {} (>= 5): polar-NURBS corner patch interior quality \
                         may be degraded; true rational Gregory for N >= 5 is pending",
                        info.n_sides
                    ),
                },
            ));
            continue;
        }
        if agg.min_achieved < requested {
            warnings.push((
                fallback_edge,
                BlendWarning::CornerContinuityDowngraded {
                    vertex: info.vertex,
                    n_sides: info.n_sides,
                    requested,
                    achieved: agg.min_achieved,
                    reason: format!(
                        "corner construction ({}) delivered {}, not {}",
                        info.method, agg.min_achieved, requested
                    ),
                },
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::CornerMethod;
    use popkern_core::topology::BRep;

    fn fabricated_patches(
        brep: &BRep,
        recipes: &[(usize, ContinuityLevel, CornerMethod)],
    ) -> Vec<CornerPatchInfo> {
        use crate::vertex_classifier::VertexBlendClassification;
        let vertex = brep.vertices.keys().next().unwrap();
        let face = brep.faces.keys().next().unwrap();
        recipes
            .iter()
            .map(|(n, cont, method)| CornerPatchInfo {
                vertex,
                n_sides: *n,
                method: *method,
                achieved_continuity: *cont,
                classification: VertexBlendClassification::SameSideConvex,
                face,
                bridge_fallbacks: Vec::new(),
                overlap_remediated: false,
                relaxed_continuity_coplanarity: None,
            })
            .collect()
    }

    #[test]
    fn corner_downgrade_fires_for_n_ge_5() {
        let brep = BRep::unit_cube(popkern_core::tolerance::Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let patches = fabricated_patches(&brep, &[(5, ContinuityLevel::G1, CornerMethod::Gregory)]);
        let opts = FilletOptions::default();
        let mut warnings = Vec::new();
        emit_corner_downgrade_warnings(&patches, &opts, edge, &mut warnings);
        assert_eq!(warnings.len(), 1, "one warning expected for N=5");
        match &warnings[0].1 {
            BlendWarning::CornerContinuityDowngraded {
                n_sides, reason, ..
            } => {
                assert_eq!(*n_sides, 5);
                assert!(
                    reason.contains("N = 5"),
                    "reason should cite N, got: {reason}"
                );
            }
            other => panic!("unexpected variant: {other:?}"),
        }
    }

    #[test]
    fn corner_downgrade_fires_for_achieved_below_requested() {
        let brep = BRep::unit_cube(popkern_core::tolerance::Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let patches = fabricated_patches(
            &brep,
            &[(3, ContinuityLevel::G1, CornerMethod::RollingBallSphere)],
        );
        let opts = FilletOptions {
            continuity: ContinuityLevel::G2,
            ..FilletOptions::default()
        };
        let mut warnings = Vec::new();
        emit_corner_downgrade_warnings(&patches, &opts, edge, &mut warnings);
        assert_eq!(warnings.len(), 1, "G1<G2 should emit one warning");
        match &warnings[0].1 {
            BlendWarning::CornerContinuityDowngraded {
                achieved,
                requested,
                ..
            } => {
                assert_eq!(*achieved, ContinuityLevel::G1);
                assert_eq!(*requested, ContinuityLevel::G2);
            }
            other => panic!("unexpected variant: {other:?}"),
        }
    }

    #[test]
    fn corner_downgrade_silent_for_matching_continuity_and_low_n() {
        let brep = BRep::unit_cube(popkern_core::tolerance::Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let patches = fabricated_patches(
            &brep,
            &[(3, ContinuityLevel::G1, CornerMethod::RollingBallSphere)],
        );
        let opts = FilletOptions::default();
        let mut warnings = Vec::new();
        emit_corner_downgrade_warnings(&patches, &opts, edge, &mut warnings);
        assert!(
            warnings.is_empty(),
            "no warning for N=3 + matching continuity"
        );
    }

    #[test]
    fn corner_downgrade_subsumes_both_conditions_into_one_warning() {
        // N>=5 AND achieved<requested — the N>=5 branch fires first
        // and emits exactly one warning; no double-report.
        let brep = BRep::unit_cube(popkern_core::tolerance::Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let patches = fabricated_patches(&brep, &[(6, ContinuityLevel::G1, CornerMethod::Gregory)]);
        let opts = FilletOptions {
            continuity: ContinuityLevel::G2,
            ..FilletOptions::default()
        };
        let mut warnings = Vec::new();
        emit_corner_downgrade_warnings(&patches, &opts, edge, &mut warnings);
        assert_eq!(warnings.len(), 1, "N>=5 subsumes the downgrade rationale");
    }

    /// R8 C7 dedup — when the three-quad dispatch produces 3
    /// `CornerPatchInfo` entries per vertex, emit exactly ONE
    /// warning, not three.
    #[test]
    fn corner_downgrade_dedups_across_three_quad_subpatches_for_same_vertex() {
        let brep = BRep::unit_cube(popkern_core::tolerance::Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        // Three patches, same vertex, same requested G2 vs achieved G1.
        let patches = fabricated_patches(
            &brep,
            &[
                (3, ContinuityLevel::G1, CornerMethod::Gregory),
                (3, ContinuityLevel::G1, CornerMethod::Gregory),
                (3, ContinuityLevel::G1, CornerMethod::Gregory),
            ],
        );
        let opts = FilletOptions {
            continuity: ContinuityLevel::G2,
            ..FilletOptions::default()
        };
        let mut warnings = Vec::new();
        emit_corner_downgrade_warnings(&patches, &opts, edge, &mut warnings);
        assert_eq!(
            warnings.len(),
            1,
            "three sub-patches at the same vertex should aggregate to one warning"
        );
    }

    /// R8 C7 dedup — when sub-patches at the same vertex report
    /// DIFFERENT achieved levels, the aggregated warning should
    /// reflect the WEAKEST (min). Synthetic coverage — today all
    /// three-quad sub-patches are G1, but future twist-extraction
    /// work may produce per-sub-patch G2 or mixed reports.
    #[test]
    fn corner_downgrade_aggregates_to_min_achieved_across_subpatches() {
        let brep = BRep::unit_cube(popkern_core::tolerance::Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let patches = fabricated_patches(
            &brep,
            &[
                (3, ContinuityLevel::G2, CornerMethod::Gregory),
                (3, ContinuityLevel::G1, CornerMethod::Gregory), // weakest
                (3, ContinuityLevel::G2, CornerMethod::Gregory),
            ],
        );
        let opts = FilletOptions {
            continuity: ContinuityLevel::G2,
            ..FilletOptions::default()
        };
        let mut warnings = Vec::new();
        emit_corner_downgrade_warnings(&patches, &opts, edge, &mut warnings);
        assert_eq!(warnings.len(), 1);
        match &warnings[0].1 {
            BlendWarning::CornerContinuityDowngraded { achieved, .. } => {
                assert_eq!(
                    *achieved,
                    ContinuityLevel::G1,
                    "aggregate should report MIN across sub-patches"
                );
            }
            other => panic!("unexpected variant: {other:?}"),
        }
    }
}
