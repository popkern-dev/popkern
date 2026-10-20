//! Interference over a subset of occurrence pairs, optionally under a
//! rigid world motion of some occurrences (PK-855): a swing checks the
//! pairs its joint moves; a film checks the parts its tracks move,
//! without touching the document's placements.

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::interference::{InterferenceOptions, PairStatus};
use crate::placement::Placement;

/// The pair's report when its base meshes PROVE an overlap (they cross
/// or nest, and `mesh_bvh::certain_overlap` holds), else `None`: the
/// uncertified check, with no distance search or refinement.
pub(crate) fn proven_overlap(
    a: &crate::interference::Occurrence,
    b: &crate::interference::Occurrence,
    tolerance: f64,
) -> Option<crate::interference::PairReport> {
    for (brep_a, mesh_a, _) in &a.bodies {
        for (brep_b, mesh_b, _) in &b.bodies {
            let meet = mesh_a.crosses(mesh_b)
                || brep_a.first_vertex().is_some_and(|p| mesh_b.contains(&p))
                || brep_b.first_vertex().is_some_and(|p| mesh_a.contains(&p));
            if meet && crate::mesh_bvh::certain_overlap(mesh_a, mesh_b, tolerance) {
                return Some(crate::interference::PairReport {
                    a: a.path.clone(),
                    b: b.path.clone(),
                    status: PairStatus::Interfering {
                        volume: None,
                        decline: None,
                    },
                    closest: None,
                    faces: None,
                });
            }
        }
    }
    None
}

/// A world-frame motion per occurrence path (`None`: in place).
pub(crate) type Motion<'a> = &'a dyn Fn(&str) -> Option<Placement>;

/// Which occurrence pairs an interference check covers, and how.
pub(crate) struct PairScope<'a> {
    /// The pairs (occurrence paths) to check.
    pub keep: &'a dyn Fn(&str, &str) -> bool,
    /// A world-frame motion applied after an occurrence's placement.
    pub motion: Option<Motion<'a>>,
    /// Refine pairs the base meshes leave open (`false`: report only
    /// the pairs the base meshes settle, a proven overlap among them).
    pub certify: bool,
}

impl PairScope<'_> {
    /// Every pair, in place, certified.
    pub(crate) fn all() -> PairScope<'static> {
        PairScope {
            keep: &|_, _| true,
            motion: None,
            certify: true,
        }
    }
}

impl Assembly {
    /// The occurrence pairs `motion` drives into each other: every pair
    /// with at least one occurrence that `moved` accepts, checked with
    /// `motion` (a world-frame motion after each occurrence's placement,
    /// `None` where it stays) applied, reported when the base meshes
    /// PROVE the overlap (`mesh_bvh::certain_overlap`). Pairs that only
    /// touch, or overlap by less than the base meshes resolve (a thread
    /// fit), are not reported: this asks whether one part was driven
    /// through another, cheaply enough to ask once per frame. The
    /// document is not changed.
    ///
    /// # Errors
    ///
    /// As [`Assembly::interference`].
    pub fn overlaps_under(
        &mut self,
        motion: &dyn Fn(&str) -> Option<Placement>,
        moved: &dyn Fn(&str) -> bool,
    ) -> Result<Vec<(String, String)>, AssemblyError> {
        let keep = |a: &str, b: &str| moved(a) || moved(b);
        let report = self.interference_where(
            &InterferenceOptions {
                overlap_volume: false,
                ..InterferenceOptions::default()
            },
            &PairScope {
                keep: &keep,
                motion: Some(motion),
                certify: false,
            },
        )?;
        Ok(report
            .pairs
            .into_iter()
            .filter(|p| matches!(p.status, PairStatus::Interfering { .. }))
            .map(|p| (p.a, p.b))
            .collect())
    }
}
