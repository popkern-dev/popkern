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

/// Append one `CornerContinuityDowngraded` warning per corner
/// patch that either (a) delivered a lower continuity than the
/// caller requested, or (b) has `n_sides >= 5` (a quality
/// caveat — single-face polar NURBS approximations degrade at
/// high valence because a degree-2 radial cross-section cannot
/// capture N-way curvature variation; true rational Gregory
/// patches for N≥5 are a §5.3 C6+ follow-up).
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
    let requested = options.continuity;
    for info in corner_patches {
        if info.n_sides >= 5 {
            warnings.push((
                fallback_edge,
                BlendWarning::CornerContinuityDowngraded {
                    vertex: info.vertex,
                    n_sides: info.n_sides,
                    requested,
                    achieved: info.achieved_continuity,
                    reason: format!(
                        "N = {} (>= 5): polar-NURBS corner patch interior quality \
                         may be degraded; true rational Gregory for N >= 5 is pending",
                        info.n_sides
                    ),
                },
            ));
            continue;
        }
        if info.achieved_continuity < requested {
            warnings.push((
                fallback_edge,
                BlendWarning::CornerContinuityDowngraded {
                    vertex: info.vertex,
                    n_sides: info.n_sides,
                    requested,
                    achieved: info.achieved_continuity,
                    reason: format!(
                        "corner construction ({}) delivered {}, not {}",
                        info.method, info.achieved_continuity, requested
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
        let vertex = brep.vertices.keys().next().unwrap();
        let face = brep.faces.keys().next().unwrap();
        recipes
            .iter()
            .map(|(n, cont, method)| CornerPatchInfo {
                vertex,
                n_sides: *n,
                method: *method,
                achieved_continuity: *cont,
                face,
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
}
