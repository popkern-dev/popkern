//! Canonical edge sampling and encroachment refinement: the stage of
//! the boundary-first tessellator that fixes every edge's samples
//! before the final face pass (split out of `tessellate.rs`).

use std::collections::{HashMap, HashSet};

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::boundary_first::EdgeSamples;
use crate::mesh::FaceMesh;
use crate::{MeshError, TessellationOptions};

/// Canonical edge samples, refined until no face's Delaunay coverage
/// asks for more, and each face's checked mesh on those samples (kept
/// for the final pass). `scope` limits both to a face subset and its
/// edges ([`crate::tessellate_faces`]); `None` is the whole body.
#[allow(clippy::type_complexity)]
pub(crate) fn sample_and_refine_edges(
    brep: &BRep,
    scope: Option<&HashSet<FaceHandle>>,
    options: &TessellationOptions,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<
    (
        HashMap<EdgeHandle, EdgeSamples>,
        HashMap<FaceHandle, Option<FaceMesh>>,
    ),
    MeshError,
> {
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(MeshError::Cancelled);
            }
        };
    }
    let cancel = tracker.token();
    // Canonical per-edge samples: adjacent faces share bit-identical
    // boundary vertices (PK-095 closure).
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let edges = scope.map(|faces| edges_of(brep, faces));
    let mut edge_samples = crate::boundary_first::sample_edges(brep, options, edges.as_ref());
    tracker.finish_stage("tessellate.edge_samples", t0, edge_samples.len());
    // Encroachment refinement (Ruppert-style): a face whose Delaunay
    // coverage fails because boundary points crowd a ring segment
    // (engraved pockets, thin annuli) requests SPLITS of the
    // canonical edge spans involved. Splits refine the SHARED sample
    // sets and all faces re-run, so watertightness is preserved by
    // construction. Bounded rounds; convergence is guaranteed by the
    // no-new-fraction check.
    const ENCROACHMENT_ROUNDS: usize = 8;
    let bf_trace = popkern_core::trace::trace_enabled("boundary-first");
    // Incremental rounds (PK-277): round 0 checks every face; later
    // rounds only faces incident to an edge whose canonical samples
    // were refined last round. A face's requests are a pure function
    // of its own boundary samples, and every skip reason below is
    // monotone (fractions only grow, per-edge caps only fill), so an
    // unchanged-boundary face can only reproduce requests that were
    // already inserted or are permanently deduped.
    let mut active: Option<std::collections::HashSet<FaceHandle>> = None;
    // Each check builds the face's boundary-first mesh, the same pure
    // function of its edge samples that the final pass would run. Keep
    // it until a split changes one of the face's edges. Rebuilding
    // every face in the final pass doubled the mesh cost of a thread
    // rod, whose rounds usually end with no requests (PK-794).
    let mut checked: HashMap<FaceHandle, Option<crate::mesh::FaceMesh>> = HashMap::new();
    for round in 0..ENCROACHMENT_ROUNDS {
        // §17 per-round poll (the per-face closures also short-circuit
        // on a fired token, so a cancelled round finishes fast).
        stage_gate!();
        let round_started = popkern_core::time::Stopwatch::start();
        // Slotmap key order, filtered — request order matches the
        // retired all-faces loop with the no-op faces removed.
        let faces: Vec<FaceHandle> = brep
            .faces
            .keys()
            .filter(|&fh| brep.faces[fh].surface.is_some())
            .filter(|fh| scope.is_none_or(|set| set.contains(fh)))
            .filter(|fh| active.as_ref().is_none_or(|set| set.contains(fh)))
            .collect();
        let checks =
            crate::face_pass::check_faces(brep, &faces, &edge_samples, options, bf_trace, cancel);
        stage_gate!();
        let mut requests = Vec::new();
        for check in checks {
            requests.extend(check.requests);
            checked.insert(check.face, check.mesh);
        }
        if requests.is_empty() {
            break;
        }
        let mut split_edges: std::collections::BTreeSet<popkern_core::topology::EdgeHandle> =
            std::collections::BTreeSet::new();
        let mut inserted = 0usize;
        for (eh, fraction) in requests {
            let Some(samples) = edge_samples.get_mut(&eh) else {
                continue;
            };
            if samples.fractions.len() >= 8192 {
                // Refinement cap per edge — matches the adaptive
                // sampler's runaway guard (96 starved the xfine
                // cut-curve splits and the coverage check bounced the
                // whole wall to the uniform-grid path; 512 refused a
                // long torus-band rim the splits its facets need,
                // PK-653).
                continue;
            }
            let idx = samples.fractions.partition_point(|&f| f < fraction);
            let near = |f: f64| (f - fraction).abs() < 1e-4;
            if (idx > 0 && near(samples.fractions[idx - 1]))
                || (idx < samples.fractions.len() && near(samples.fractions[idx]))
            {
                continue; // already refined here
            }
            let Some(curve) = brep.edges.get(eh).and_then(|e| e.curve.as_ref()) else {
                continue;
            };
            use popkern_core::geometry::curve::CurveEval;
            let ceval: &dyn CurveEval = &curve.geometry;
            let point = ceval.point_at(curve.t_start + (curve.t_end - curve.t_start) * fraction);
            samples.fractions.insert(idx, fraction);
            samples.points.insert(idx, point);
            split_edges.insert(eh);
            inserted += 1;
        }
        if bf_trace && inserted > 0 {
            eprintln!("  [bf] encroachment round {round}: split {inserted} edge spans");
        }
        tracker.finish_stage("tessellate.encroachment_round", round_started, inserted);
        if inserted == 0 {
            break;
        }
        let mut changed = crate::face_pass::faces_incident_to_edges(brep, &split_edges);
        if let Some(scope) = scope {
            changed.retain(|f| scope.contains(f));
        }
        for face in &changed {
            checked.remove(face);
        }
        active = Some(changed);
    }
    Ok((edge_samples, checked))
}

/// Every edge bounding one of `faces`.
fn edges_of(brep: &BRep, faces: &HashSet<FaceHandle>) -> HashSet<EdgeHandle> {
    let mut edges = HashSet::new();
    for &fh in faces {
        let Some(face) = brep.faces.get(fh) else {
            continue;
        };
        for &lh in std::iter::once(&face.outer_loop).chain(&face.inner_loops) {
            let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
                continue;
            };
            let mut he = start;
            // Bounded walk: a malformed ring must not spin forever.
            for _ in 0..brep.half_edges.len() {
                let Some(h) = brep.half_edges.get(he) else {
                    break;
                };
                edges.insert(h.edge);
                he = h.next;
                if he == start {
                    break;
                }
            }
        }
    }
    edges
}

/// Meshes of `faces` of `brep` alone, by the whole-body pipeline's
/// boundary-first face pass: canonical samples on the faces' own edges,
/// encroachment refined among those faces, the same per-face fallback.
/// It is a measuring mesh, not a display one: the cross-face finish
/// (projection sewing, coherent orientation, facet-tilt splits) is not
/// run, so callers bound geometry by each triangle's own deviation.
/// [`crate::tessellate_face`] alone meshes a face on its uniform
/// subdivision grid, ~500 times the triangles at a fine chord (a
/// certification region, PK-900).
///
/// # Errors
///
/// The first face's tessellation error, in `faces` order.
pub fn tessellate_faces(
    brep: &BRep,
    faces: &[FaceHandle],
    options: &TessellationOptions,
) -> Result<Vec<FaceMesh>, MeshError> {
    let scope: HashSet<FaceHandle> = faces.iter().copied().collect();
    let mut tracker = popkern_core::progress::StageTracker::inert();
    let (edge_samples, checked) =
        sample_and_refine_edges(brep, Some(&scope), options, &mut tracker)?;
    let faces: Vec<FaceHandle> = faces
        .iter()
        .copied()
        .filter(|&fh| brep.faces.get(fh).is_some_and(|f| f.surface.is_some()))
        .collect();
    let trace = popkern_core::trace::trace_enabled("boundary-first");
    crate::face_pass::tessellate_face_pass(
        brep,
        &faces,
        checked,
        &edge_samples,
        options,
        trace,
        None,
    )
    .into_iter()
    .map(|outcome| match outcome {
        crate::face_pass::FacePassOutcome::BoundaryFirst(fm) => Ok(fm),
        crate::face_pass::FacePassOutcome::UniformGrid { mesh, .. } => mesh,
    })
    .collect()
}

#[cfg(test)]
#[path = "encroachment_tests.rs"]
mod tests;
