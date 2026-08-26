//! Skin collapse detection: decline offsets that fold the skin.
//!
//! After vertex/edge recompute the offset skin's vertices are final. If
//! the requested distance exceeds the local feature size, edges shrink
//! through zero and invert — the skin folds through itself and the
//! result is inside-out wrong geometry that can still be structurally
//! consistent (wrong-shape-served-clean class). This check compares
//! every recomputed edge chord against its original: a chord whose
//! projection onto the original direction is not positive has collapsed,
//! and the operation declines with [`ShellError::OffsetTooLarge`] naming
//! the binding edge and an estimated maximum safe distance.
//!
//! The estimate uses the linear model of vertex motion in the offset
//! distance (exact for plane triple-point vertices, first-order
//! otherwise): a chord shrinking from `len0` to `proj` at distance `d`
//! reaches zero at `d* = d·len0/(len0 − proj)`; the reported maximum is
//! `0.99·d*`.
//!
//! Closed edges (chord below tolerance in the ORIGINAL body) carry no
//! usable chord direction and are skipped; curved-skin folding beyond
//! chord inversion is the trim stage's job
//! (docs/architecture/solid-offset-shell-design.md §5).

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// Routing verdict of the chord gate (PK-327): tangential zero-length
/// collapses decline immediately, but INVERTED chords (the skin edge
/// crossed through zero — a transversal collision) are routed to the
/// global self-intersection trim stage with the decline deferred. If
/// the trim stage cannot serve, the deferred error is returned so
/// behavior degrades to the honest pre-trim decline.
#[derive(Debug, Default)]
pub struct TrimRouting {
    /// Outer-body edges whose skin image inverted (crossed zero).
    pub inverted_edges: Vec<EdgeHandle>,
    /// The decline the chord gate WOULD have issued, deferred while
    /// the trim stage decides.
    pub deferred: Option<ShellError>,
}

/// Decline when any recomputed skin edge has collapsed or inverted.
///
/// `outer` is the original body, `skin` the recomputed offset skin, and
/// `skin_map` the outer→skin handle map. `requested` is the caller's
/// signed offset/thickness parameter (used for reporting and the linear
/// max-safe estimate).
///
/// # Errors
///
/// [`ShellError::OffsetTooLarge`] on the first TANGENTIAL collapse (an
/// edge shrinking to zero without crossing — a degenerate contact the
/// trim stage cannot resolve). Inverted edges are returned in the
/// routing instead, with the matching decline deferred.
pub fn check_skin_collapse(
    outer: &BRep,
    skin: &BRep,
    skin_map: &SkinMap,
    requested: f64,
    tolerance: &Tolerance,
) -> Result<TrimRouting, ShellError> {
    let mut routing = TrimRouting::default();
    let tol = tolerance.linear();
    // Deterministic order: iterate the outer body's edge slotmap (stable
    // insertion order), not the HashMap.
    for (outer_eh, edge) in &outer.edges {
        let Some(&skin_eh) = skin_map.edges.get(&outer_eh) else {
            continue; // edge of a removed face
        };
        if !skin.edges.contains_key(skin_eh) {
            continue;
        }

        let he = &outer.half_edges[edge.half_edge];
        let twin = &outer.half_edges[he.twin];
        let p0 = outer.vertices[he.origin].position;
        let p1 = outer.vertices[twin.origin].position;
        let old_chord = p1 - p0;
        let len0 = old_chord.length();
        if len0 <= tol {
            continue; // closed edge: no chord direction to compare
        }
        let old_dir = old_chord * (1.0 / len0);

        let skin_edge = &skin.edges[skin_eh];
        let she = &skin.half_edges[skin_edge.half_edge];
        let stwin = &skin.half_edges[she.twin];
        let q0 = skin.vertices[she.origin].position;
        let q1 = skin.vertices[stwin.origin].position;
        let proj = (q1 - q0).dot(&old_dir);

        if proj <= tol {
            let magnitude = requested.abs();
            // Linear vertex-motion model: chord length hits zero at d*.
            let shrink = len0 - proj;
            let max_safe = if shrink > tol {
                0.99 * magnitude * len0 / shrink
            } else {
                0.0
            };
            let decline = ShellError::OffsetTooLarge {
                requested: magnitude,
                max_safe,
                binding: format!("edge {outer_eh:?} (length {len0:.6}) collapses at this distance"),
            };
            if proj < -tol {
                // Inverted: the edge crossed through zero — a
                // transversal collision the trim stage can remove.
                routing.inverted_edges.push(outer_eh);
                if routing.deferred.is_none() {
                    routing.deferred = Some(decline);
                }
            } else {
                // Tangential zero-length contact: degenerate, decline.
                return Err(decline);
            }
        }
    }
    Ok(routing)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::classify_faces;
    use crate::inner_skin::build_inner_skin;
    use crate::surface_offset::compute_offset_surfaces;
    use crate::types::ShellDirection;
    use crate::vertex_recompute::recompute_vertices;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    fn recomputed_inner(thickness: f64) -> (BRep, BRep, SkinMap) {
        let outer = BRep::unit_cube(tol()).unwrap();
        let classification =
            classify_faces(&outer, thickness, &[], ShellDirection::Inward, &[]).unwrap();
        let offset = compute_offset_surfaces(&outer, &classification, true, &tol()).unwrap();
        let (mut inner, skin_map) = build_inner_skin(&outer, &[]).unwrap();
        recompute_vertices(&mut inner, &offset.surfaces, &skin_map, &outer, &tol()).unwrap();
        (outer, inner, skin_map)
    }

    #[test]
    fn thin_shell_passes() {
        let (outer, inner, map) = recomputed_inner(0.1);
        check_skin_collapse(&outer, &inner, &map, 0.1, &tol()).unwrap();
    }

    #[test]
    fn crossing_shell_routes_inverted_edges_with_deferred_decline() {
        // Unit cube at 0.6 inward: opposite faces cross, every edge
        // inverts — routed to the trim stage with the decline DEFERRED
        // (PK-327), not returned.
        let (outer, inner, map) = recomputed_inner(0.6);
        let routing = check_skin_collapse(&outer, &inner, &map, 0.6, &tol()).unwrap();
        assert!(!routing.inverted_edges.is_empty(), "inverted edges routed");
        match routing.deferred {
            Some(ShellError::OffsetTooLarge {
                requested,
                max_safe,
                binding,
            }) => {
                assert!((requested - 0.6).abs() < 1e-12);
                assert!(
                    max_safe > 0.0 && max_safe < 0.6,
                    "max_safe {max_safe} must sit below the requested 0.6"
                );
                assert!(binding.contains("edge"), "{binding}");
            }
            other => panic!("expected deferred OffsetTooLarge, got {other:?}"),
        }
    }

    #[test]
    fn exactly_degenerate_shell_declines() {
        // At 0.5 the inner cube collapses to a point: proj == 0.
        let (outer, inner, map) = recomputed_inner(0.5);
        assert!(check_skin_collapse(&outer, &inner, &map, 0.5, &tol()).is_err());
    }

    #[test]
    fn max_safe_estimate_is_serveable() {
        // The deferred max_safe from a routed 0.6 shell must itself pass.
        let (outer, _, _) = recomputed_inner(0.6);
        let (_, inner, map) = recomputed_inner(0.6);
        let routing = check_skin_collapse(&outer, &inner, &map, 0.6, &tol()).unwrap();
        let Some(ShellError::OffsetTooLarge { max_safe, .. }) = routing.deferred else {
            panic!("expected deferred OffsetTooLarge");
        };
        let routing2 = {
            let (outer2, inner2, map2) = recomputed_inner(max_safe);
            check_skin_collapse(&outer2, &inner2, &map2, max_safe, &tol()).unwrap()
        };
        assert!(routing2.deferred.is_none(), "max_safe must serve clean");
        assert!(routing2.inverted_edges.is_empty());
    }
}
