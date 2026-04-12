//! Deferred ON-face resolution for ambiguous classifications (§6.4 [BEYOND]).
//!
//! Today the boolean engine resolves ON classifications immediately:
//! [`crate::classify_on_face::classify_on_face`] only processes face
//! overlaps with `full_overlap == true`, and partial overlaps fall
//! through to ray-casting. Ray-casting on a partial overlap is
//! tolerance-sensitive and order-dependent, which is exactly the
//! cliff edge §6.4 calls out:
//!
//! > "When ON classification is ambiguous (tolerance-induced or
//! > partial overlap), classify the face as ON and defer resolution
//! > to a topology-aware disambiguation phase that considers the
//! > global assembly context — neighboring faces, shell closure
//! > requirements, and orientation consistency."
//!
//! This module supplies the infrastructure for that disambiguation
//! phase. It does not yet replace the immediate-resolution path —
//! callers opt in by running [`find_ambiguous_pairs`] before the
//! main classification, and then [`disambiguate_with_neighbors`]
//! after ray-casting has produced its initial answer.
//!
//! The disambiguation rule today is local: a partial-overlap face
//! whose neighbors classify as `Inside` should also be `Inside`,
//! and likewise for `Outside`. Future work can extend this to
//! shell-closure-aware reasoning.

use popkern_core::topology::{BRep, FaceHandle};

use crate::coincidence::types::{CoincidenceMap, CoincidenceOrientation, FaceOverlap};
use crate::error::Classification;

/// A face pair that classify_on_face cannot resolve immediately.
///
/// `face_a` and `face_b` are coincident in the partial-overlap
/// sense (`full_overlap == false`); the orientation is preserved
/// from the source overlap so the disambiguation phase can apply
/// the orientation-aware rules.
#[derive(Debug, Clone)]
pub struct AmbiguousOnPair {
    /// Face from operand A.
    pub face_a: FaceHandle,
    /// Face from operand B.
    pub face_b: FaceHandle,
    /// Orientation relationship at the overlap region.
    pub orientation: CoincidenceOrientation,
}

/// Walk the coincidence map and collect every face pair that
/// classify_on_face would skip — i.e. partial overlaps.
///
/// The returned vector is in the order the overlaps appear in
/// `coincidence.face_overlaps`, so the disambiguation phase can
/// match it against the source overlap if needed.
pub fn find_ambiguous_pairs(coincidence: &CoincidenceMap) -> Vec<AmbiguousOnPair> {
    coincidence
        .face_overlaps
        .iter()
        .filter(|fo| !fo.full_overlap)
        .map(|fo: &FaceOverlap| AmbiguousOnPair {
            face_a: fo.face_a,
            face_b: fo.face_b,
            orientation: fo.orientation,
        })
        .collect()
}

/// Audit a classification result for ambiguous-ON disagreements.
///
/// Walks every ambiguous pair from `coincidence` and, for each one
/// that has a corresponding entry in `face_classifications`, runs
/// [`disambiguate_with_neighbors`] using the other classifications
/// as neighbour votes. Returns a list of `(face_handle, original,
/// suggested)` tuples where the disambiguation would change the
/// classification — useful as a diagnostic to see how often the
/// deferred resolver would override the immediate decision.
///
/// The audit is read-only — it never mutates the classifications.
/// Callers can either ignore the report or apply it themselves.
pub fn audit_ambiguous_on_classifications(
    brep: &BRep,
    coincidence: &CoincidenceMap,
    face_classifications: &[(FaceHandle, Classification)],
) -> Vec<(FaceHandle, Classification, Classification)> {
    let pairs = find_ambiguous_pairs(coincidence);
    let mut disagreements = Vec::new();
    for pair in &pairs {
        // Try both face_a and face_b in case the classification map
        // only contains one side of the pair.
        for face in [pair.face_a, pair.face_b] {
            if let Some((_, original)) = face_classifications.iter().find(|(f, _)| *f == face) {
                let suggested =
                    disambiguate_with_neighbors(face, *original, brep, face_classifications);
                if suggested != *original {
                    disagreements.push((face, *original, suggested));
                }
            }
        }
    }
    disagreements
}

/// Resolve an ambiguous-ON classification using neighbour context.
///
/// Given the ray-casting result for `face` and the classifications
/// of every adjacent face in the same shell, returns a resolved
/// classification that is consistent with the dominant neighbour
/// class. If `face` has no neighbours, returns the ray-casting
/// result unchanged.
///
/// The disambiguation rule:
/// - If a strict majority of neighbours are `Inside`, the face is
///   `Inside`.
/// - If a strict majority are `Outside`, the face is `Outside`.
/// - Otherwise, the original `ray_classification` survives.
///
/// This is a local rule, not the full topology-aware reasoning §6.4
/// [BEYOND] envisions; it eliminates the most common
/// near-tolerance flips while leaving harder cases for future work.
pub fn disambiguate_with_neighbors(
    face: FaceHandle,
    ray_classification: Classification,
    brep: &BRep,
    neighbor_classifications: &[(FaceHandle, Classification)],
) -> Classification {
    let neighbors = collect_face_neighbors(brep, face);
    if neighbors.is_empty() {
        return ray_classification;
    }

    let mut inside_votes = 0_usize;
    let mut outside_votes = 0_usize;
    for &(nf, class) in neighbor_classifications {
        if !neighbors.contains(&nf) {
            continue;
        }
        match class {
            Classification::Inside => inside_votes += 1,
            Classification::Outside => outside_votes += 1,
            // ON-classified neighbors carry no inside/outside vote;
            // they are themselves ambiguous.
            Classification::OnBoundaryCooriented | Classification::OnBoundaryCrossing => {}
        }
    }

    let total = inside_votes + outside_votes;
    if total == 0 {
        return ray_classification;
    }
    if inside_votes * 2 > total {
        Classification::Inside
    } else if outside_votes * 2 > total {
        Classification::Outside
    } else {
        ray_classification
    }
}

/// Collect the face handles directly adjacent to `face` via shared
/// edges in `brep`. Used by [`disambiguate_with_neighbors`].
fn collect_face_neighbors(brep: &BRep, face: FaceHandle) -> Vec<FaceHandle> {
    let mut neighbors = Vec::new();
    let Some(face_data) = brep.faces.get(face) else {
        return neighbors;
    };
    let outer_loop = face_data.outer_loop;
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return neighbors;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        let twin = he.twin;
        if let Some(twin_he) = brep.half_edges.get(twin)
            && let Some(twin_loop) = brep.loops.get(twin_he.loop_)
        {
            let nf = twin_loop.face;
            if nf != face && !neighbors.contains(&nf) {
                neighbors.push(nf);
            }
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    neighbors
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coincidence::types::{CoincidenceOrientation, FaceOverlap};
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    fn make_overlap(full: bool) -> FaceOverlap {
        FaceOverlap {
            face_a: FaceHandle::default(),
            face_b: FaceHandle::default(),
            orientation: CoincidenceOrientation::Same,
            full_overlap: full,
            overlap_center: Point3::new(0.0, 0.0, 0.0),
        }
    }

    #[test]
    fn find_ambiguous_pairs_skips_full_overlaps() {
        let mut map = CoincidenceMap::default();
        map.face_overlaps.push(make_overlap(true));
        map.face_overlaps.push(make_overlap(false));
        map.face_overlaps.push(make_overlap(true));
        let pairs = find_ambiguous_pairs(&map);
        assert_eq!(
            pairs.len(),
            1,
            "only the partial-overlap pair should be ambiguous"
        );
    }

    #[test]
    fn find_ambiguous_pairs_empty_map() {
        let map = CoincidenceMap::default();
        let pairs = find_ambiguous_pairs(&map);
        assert!(pairs.is_empty());
    }

    #[test]
    fn disambiguate_with_no_neighbors_returns_ray_class() {
        let cube = BRep::unit_cube(tol()).unwrap();
        let face = cube.faces.keys().next().unwrap();
        // Empty neighbor classifications → no votes → ray result survives.
        let resolved = disambiguate_with_neighbors(face, Classification::Inside, &cube, &[]);
        assert_eq!(resolved, Classification::Inside);
    }

    #[test]
    fn disambiguate_inside_majority_picks_inside() {
        let cube = BRep::unit_cube(tol()).unwrap();
        let face = cube.faces.keys().next().unwrap();
        let neighbors = collect_face_neighbors(&cube, face);
        // Vote every neighbour as Inside; the face was originally Outside;
        // the majority rule should override to Inside.
        let votes: Vec<(FaceHandle, Classification)> = neighbors
            .iter()
            .map(|&n| (n, Classification::Inside))
            .collect();
        let resolved = disambiguate_with_neighbors(face, Classification::Outside, &cube, &votes);
        assert_eq!(resolved, Classification::Inside);
    }

    #[test]
    fn audit_with_no_ambiguous_pairs_returns_empty() {
        let cube = BRep::unit_cube(tol()).unwrap();
        let face = cube.faces.keys().next().unwrap();
        let coincidence = CoincidenceMap::default();
        let classes = vec![(face, Classification::Inside)];
        let report = audit_ambiguous_on_classifications(&cube, &coincidence, &classes);
        assert!(report.is_empty());
    }

    #[test]
    fn audit_returns_disagreement_when_neighbors_outvote() {
        let cube = BRep::unit_cube(tol()).unwrap();
        let target = cube.faces.keys().next().unwrap();

        // Vote every neighbor as Inside; the target is Outside.
        // The audit should disagree with the original classification.
        let neighbors = collect_face_neighbors(&cube, target);
        let mut classes: Vec<(FaceHandle, Classification)> = neighbors
            .iter()
            .map(|&n| (n, Classification::Inside))
            .collect();
        classes.push((target, Classification::Outside));

        // Build a coincidence map with target as a partial overlap so
        // it's flagged as ambiguous.
        let mut coincidence = CoincidenceMap::default();
        coincidence.face_overlaps.push(make_overlap(false));
        // The default overlap uses FaceHandle::default(); fix the
        // face_a to be our actual target so the audit picks it up.
        coincidence.face_overlaps[0].face_a = target;

        let report = audit_ambiguous_on_classifications(&cube, &coincidence, &classes);
        assert!(
            !report.is_empty(),
            "audit should report disagreement when neighbors outvote target"
        );
        let (face, original, suggested) = report[0];
        assert_eq!(face, target);
        assert_eq!(original, Classification::Outside);
        assert_eq!(suggested, Classification::Inside);
    }

    #[test]
    fn disambiguate_split_vote_keeps_ray_class() {
        let cube = BRep::unit_cube(tol()).unwrap();
        let face = cube.faces.keys().next().unwrap();
        let neighbors = collect_face_neighbors(&cube, face);
        // Half inside, half outside → no strict majority → ray result wins.
        let mut votes: Vec<(FaceHandle, Classification)> = Vec::new();
        for (i, &n) in neighbors.iter().enumerate() {
            votes.push((
                n,
                if i % 2 == 0 {
                    Classification::Inside
                } else {
                    Classification::Outside
                },
            ));
        }
        let resolved = disambiguate_with_neighbors(face, Classification::Outside, &cube, &votes);
        assert_eq!(resolved, Classification::Outside);
    }
}
