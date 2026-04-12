//! Lightweight interference detection between two solids.
//!
//! Implements §2.4 of the boolean requirements: a fast "do these two
//! bodies overlap?" check that runs in a small fraction of the cost of
//! a full boolean. Reuses the BVH-accelerated AABB filter from the
//! intersection phase, then runs surface-surface intersection on
//! candidate pairs only until enough evidence is gathered.
//!
//! The check has two modes:
//!
//! - [`InterferenceMode::FirstHit`] terminates as soon as the first
//!   genuinely intersecting face pair is confirmed. Use this for
//!   "do these touch?" yes/no queries.
//! - [`InterferenceMode::AllPairs`] visits every BVH-candidate pair
//!   and returns the full list of interfering face pairs. Use this for
//!   visualization, clash highlighting, or as a diagnostic when
//!   first-hit isn't enough.
//!
//! Both modes return an interference AABB (the union of all measured
//! overlap regions) and a `pairs_tested` counter for telemetry.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::intersection::{SurfaceSurfaceIntersection, intersect_surfaces};
use popkern_core::geometry::uv_projection::bounded_parameter_range;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::face_vertex_positions;
use crate::intersect::compute_face_aabbs;

/// Mode for [`interference_check`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterferenceMode {
    /// Terminate at the first confirmed face-pair intersection.
    /// Returns immediately once `interfering` is known to be true.
    FirstHit,
    /// Visit every BVH-candidate pair and return all interfering pairs.
    AllPairs,
}

/// Report from an interference check.
///
/// Always populated, even when no interference is found. Inspect
/// [`Self::interfering`] for the yes/no answer; the other fields
/// describe the search effort and (when interference exists) where
/// it occurs.
#[derive(Debug, Clone, Default)]
pub struct InterferenceReport {
    /// True iff at least one face pair was confirmed to intersect.
    pub interfering: bool,
    /// Pairs whose surfaces were confirmed to intersect.
    /// In `FirstHit` mode contains exactly one entry; in `AllPairs`
    /// mode contains every confirmed pair.
    pub interfering_pairs: Vec<(FaceHandle, FaceHandle)>,
    /// Union of overlap AABBs across all confirmed pairs. `None` if
    /// no interference was found.
    pub interference_aabb: Option<Aabb>,
    /// Number of candidate face pairs the BVH returned for testing.
    pub candidate_pairs: usize,
    /// Number of candidate pairs actually tested via surface-surface
    /// intersection. `pairs_tested <= candidate_pairs`; in
    /// `FirstHit` mode this can be much smaller because the loop
    /// terminates early.
    pub pairs_tested: usize,
}

/// Run an interference check between two BRep solids.
///
/// Returns an [`InterferenceReport`] describing whether the two
/// bodies overlap, where, and how much work was done. See
/// [`InterferenceMode`] for the cost/completeness tradeoff.
///
/// # Performance
///
/// `FirstHit` mode is sub-linear in the number of face pairs for
/// disjoint bodies (BVH rejects everything) and bounded by O(k) for
/// touching bodies, where `k` is the number of pairs that share an
/// AABB. Worst case it equals `AllPairs`. `AllPairs` is bounded by
/// the number of BVH-candidate pairs, which is much smaller than the
/// quadratic `n_a * n_b` for typical models.
pub fn interference_check(
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
    mode: InterferenceMode,
) -> InterferenceReport {
    let aabbs_a = compute_face_aabbs(a);
    let aabbs_b = compute_face_aabbs(b);
    let n_a = aabbs_a.len();

    // Build a BVH over both bodies' face AABBs (matching `intersect.rs`).
    use popkern_core::validation::bvh::Bvh;
    let mut combined: Vec<(usize, Aabb)> = Vec::with_capacity(n_a + aabbs_b.len());
    for (i, (_, aabb)) in aabbs_a.iter().enumerate() {
        combined.push((i, aabb.clone()));
    }
    for (i, (_, aabb)) in aabbs_b.iter().enumerate() {
        combined.push((n_a + i, aabb.clone()));
    }
    let bvh = Bvh::build(&combined);
    let candidate_pairs = bvh.find_overlapping_pairs(tolerance.linear());

    let mut report = InterferenceReport {
        candidate_pairs: candidate_pairs.len(),
        ..InterferenceReport::default()
    };

    for (id_x, id_y) in &candidate_pairs {
        // BVH may emit pairs from the same body — only cross-body pairs
        // matter for interference between A and B.
        let (a_idx, b_idx) = if *id_x < n_a && *id_y >= n_a {
            (*id_x, *id_y - n_a)
        } else if *id_y < n_a && *id_x >= n_a {
            (*id_y, *id_x - n_a)
        } else {
            continue;
        };

        let (face_a, aabb_a) = &aabbs_a[a_idx];
        let (face_b, aabb_b) = &aabbs_b[b_idx];

        let surf_a = match a.faces[*face_a].surface.as_ref() {
            Some(s) => &s.geometry,
            None => continue,
        };
        let surf_b = match b.faces[*face_b].surface.as_ref() {
            Some(s) => &s.geometry,
            None => continue,
        };

        let pos_a = face_vertex_positions(a, *face_a);
        let pos_b = face_vertex_positions(b, *face_b);
        let (u1, v1) = bounded_parameter_range(&pos_a, surf_a);
        let (u2, v2) = bounded_parameter_range(&pos_b, surf_b);

        report.pairs_tested += 1;

        let result = match intersect_surfaces(surf_a, u1, v1, surf_b, u2, v2, tolerance) {
            Ok(r) => r,
            Err(_) => continue,
        };

        // A non-empty intersection result means the surfaces actually
        // touch within their bounded UV ranges. Both Curve and
        // CoincidentRegion variants count as interference.
        let hit = result.intersections.iter().any(|isect| {
            matches!(
                isect,
                SurfaceSurfaceIntersection::Curve(_)
                    | SurfaceSurfaceIntersection::CoincidentRegion { .. }
            )
        });

        if hit {
            report.interfering = true;
            report.interfering_pairs.push((*face_a, *face_b));
            // Union of (face_a aabb ∩ face_b aabb) into the report
            // bbox. This is a rough overlap region — accurate enough
            // for visualization and clash-highlighting, much cheaper
            // than computing actual intersection-curve bounding boxes.
            let overlap = intersect_aabbs(aabb_a, aabb_b);
            report.interference_aabb = Some(match report.interference_aabb.take() {
                Some(existing) => union_aabbs(&existing, &overlap),
                None => overlap,
            });

            if mode == InterferenceMode::FirstHit {
                return report;
            }
        }
    }

    report
}

/// Compute the AABB of the intersection of two AABBs.
///
/// Caller must ensure the inputs overlap; otherwise the returned AABB
/// is degenerate (`max < min` along at least one axis).
fn intersect_aabbs(a: &Aabb, b: &Aabb) -> Aabb {
    use popkern_core::geometry::Point3;
    Aabb {
        min: Point3::new(
            a.min.x.max(b.min.x),
            a.min.y.max(b.min.y),
            a.min.z.max(b.min.z),
        ),
        max: Point3::new(
            a.max.x.min(b.max.x),
            a.max.y.min(b.max.y),
            a.max.z.min(b.max.z),
        ),
    }
}

/// Compute the AABB enclosing both inputs.
fn union_aabbs(a: &Aabb, b: &Aabb) -> Aabb {
    use popkern_core::geometry::Point3;
    Aabb {
        min: Point3::new(
            a.min.x.min(b.min.x),
            a.min.y.min(b.min.y),
            a.min.z.min(b.min.z),
        ),
        max: Point3::new(
            a.max.x.max(b.max.x),
            a.max.y.max(b.max.y),
            a.max.z.max(b.max.z),
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Transform3;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn disjoint_cubes_report_no_interference() {
        let a = BRep::unit_cube(tol()).unwrap();
        let b = BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(5.0, 0.0, 0.0))
            .unwrap();
        let report = interference_check(&a, &b, &tol(), InterferenceMode::FirstHit);
        assert!(!report.interfering);
        assert!(report.interfering_pairs.is_empty());
        assert!(report.interference_aabb.is_none());
    }

    #[test]
    fn overlapping_cubes_report_interference() {
        let a = BRep::unit_cube(tol()).unwrap();
        let b = BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(0.5, 0.3, 0.2))
            .unwrap();
        let report = interference_check(&a, &b, &tol(), InterferenceMode::FirstHit);
        assert!(report.interfering);
        assert!(!report.interfering_pairs.is_empty());
        assert!(report.interference_aabb.is_some());
    }

    #[test]
    fn first_hit_visits_fewer_pairs_than_all() {
        let a = BRep::unit_cube(tol()).unwrap();
        let b = BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(0.5, 0.3, 0.2))
            .unwrap();
        let first = interference_check(&a, &b, &tol(), InterferenceMode::FirstHit);
        let all = interference_check(&a, &b, &tol(), InterferenceMode::AllPairs);
        assert!(first.interfering && all.interfering);
        // FirstHit visits at most as many pairs as AllPairs.
        assert!(first.pairs_tested <= all.pairs_tested);
        // For overlapping cubes, AllPairs finds multiple interfering pairs.
        assert!(all.interfering_pairs.len() >= first.interfering_pairs.len());
    }

    #[test]
    fn touching_cubes_report_interference() {
        // Two unit cubes sharing a single face.
        let a = BRep::unit_cube(tol()).unwrap();
        let b = BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(1.0, 0.0, 0.0))
            .unwrap();
        let report = interference_check(&a, &b, &tol(), InterferenceMode::AllPairs);
        // Touching cubes should be detected as interfering (within tolerance).
        assert!(
            report.interfering,
            "touching cubes should report interference, pairs_tested={}, candidates={}",
            report.pairs_tested, report.candidate_pairs
        );
    }
}
