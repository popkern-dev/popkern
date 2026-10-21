//! Per-face tessellation pass dispatch (PK-277).
//!
//! Both per-face passes of `tessellate_brep` — encroachment
//! split-request collection and final face meshing — are pure
//! functions of `&BRep` (plain `SlotMap` storage, no interior
//! mutability), the shared read-only canonical edge samples, and the
//! options, so faces can run on separate threads. Per-face results are
//! collected positionally (rayon's indexed `collect`) and flattened in
//! caller face order, so output is identical at any thread count.
//! the `boundary-first` trace topic forces the sequential path so interleaved thread
//! output cannot garble the dump.

use std::collections::{BTreeSet, HashMap, HashSet};

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::boundary_first::{EdgeSamples, SplitRequest};
use crate::mesh::FaceMesh;
use crate::{MeshError, TessellationOptions};

/// Outcome of the final tessellation pass for one face.
pub(crate) enum FacePassOutcome {
    /// The boundary-first path served the face.
    BoundaryFirst(FaceMesh),
    /// The face fell back to the uniform-grid pcurve path.
    UniformGrid {
        /// True when the face's outer loop is a pure seam spine
        /// (every edge traversed twice) — such faces keep their exact
        /// chart sampling and are excluded from projection sewing.
        seam_spine: bool,
        /// The uniform-grid tessellation result; the caller propagates the
        /// first error in face order.
        mesh: Result<FaceMesh, MeshError>,
    },
}

/// One face's encroachment check: the split requests it raised and
/// the boundary-first mesh built alongside them (`None` when the
/// boundary-first path declined and the face needs the uniform-grid
/// fallback).
pub(crate) struct FaceCheck {
    /// The face checked.
    pub(crate) face: FaceHandle,
    /// Its split requests, in the order the face raised them.
    pub(crate) requests: Vec<SplitRequest>,
    /// The boundary-first mesh, before the final pass's normal repair
    /// and angular refinement.
    pub(crate) mesh: Option<FaceMesh>,
}

/// Run the encroachment check for `faces`, in face order.
///
/// Each face is checked independently, so flattening the requests in
/// face order gives the same sequence as running the faces one by one.
/// The caller's insert-and-dedup pass sees the same requests in the
/// same order at any thread count. A face skipped by a fired cancel
/// token is left out; the caller's round poll then discards the run.
pub(crate) fn check_faces(
    brep: &BRep,
    faces: &[FaceHandle],
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
    options: &TessellationOptions,
    trace: bool,
    cancel: Option<&popkern_core::cancellation::CancellationToken>,
) -> Vec<FaceCheck> {
    let per_face = |face_handle: FaceHandle| -> Option<FaceCheck> {
        // PK-276: a fired token short-circuits the remaining per-face
        // work; the caller's round poll surfaces the cancellation, and
        // a cancelled run's output is discarded wholesale, so the
        // truncated request list never reaches a served mesh.
        if cancel.is_some_and(|t| t.is_cancelled()) {
            return None;
        }
        let mut requests = Vec::new();
        let mesh = crate::boundary_first::tessellate_face_boundary_first(
            brep,
            face_handle,
            edge_samples,
            options,
            &mut requests,
        );
        Some(FaceCheck {
            face: face_handle,
            requests,
            mesh,
        })
    };
    dispatch(faces.to_vec(), trace, per_face)
        .into_iter()
        .flatten()
        .collect()
}

/// Run the final tessellation pass for `faces`, in face order.
///
/// Boundary-first is attempted per face; fallbacks run the uniform-grid path
/// with canonical samples. Uniform-grid errors are carried in the outcome so
/// the caller can propagate the first one in face order — positionally
/// identical to the retired sequential loop's early `?` return.
///
/// `checked` holds the boundary-first result of every face whose
/// encroachment check ran on the final edge samples. Such a face is not
/// rebuilt: its boundary-first mesh is a pure function of its own
/// edges' samples, which have not changed since the check (PK-794).
pub(crate) fn tessellate_face_pass(
    brep: &BRep,
    faces: &[FaceHandle],
    mut checked: HashMap<FaceHandle, Option<FaceMesh>>,
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
    options: &TessellationOptions,
    trace: bool,
    cancel: Option<&popkern_core::cancellation::CancellationToken>,
) -> Vec<FacePassOutcome> {
    let items: Vec<(FaceHandle, Option<Option<FaceMesh>>)> = faces
        .iter()
        .map(|&face_handle| (face_handle, checked.remove(&face_handle)))
        .collect();
    let per_face =
        |(face_handle, checked): (FaceHandle, Option<Option<FaceMesh>>)| -> FacePassOutcome {
            // PK-276: a fired token short-circuits remaining faces with a
            // Cancelled outcome; the caller polls before consuming, so the
            // whole pass surfaces as MeshError::Cancelled.
            if cancel.is_some_and(|t| t.is_cancelled()) {
                return FacePassOutcome::UniformGrid {
                    seam_spine: false,
                    mesh: Err(MeshError::Cancelled),
                };
            }
            let boundary_first = checked.unwrap_or_else(|| {
                crate::boundary_first::tessellate_face_boundary_first(
                    brep,
                    face_handle,
                    edge_samples,
                    options,
                    &mut Vec::new(),
                )
            });
            if let Some(fm) = boundary_first {
                let mut fm = fm;
                if let Some(sof) = brep.faces.get(face_handle).and_then(|f| f.surface.as_ref()) {
                    crate::vertex_normal::repair_singular_normals(&mut fm, &sof.geometry);
                    // Angular facet bound: boundary-first meshes only (the
                    // uniform-grid fallback below stays chord-only — see
                    // `tessellate_face_with_samples`).
                    if let Some(max_angle) = options.angular_bound() {
                        crate::angular_refine::refine_face_mesh(&mut fm, &sof.geometry, max_angle);
                        crate::tilt_refine::finish_tilt(
                            &mut fm,
                            &sof.geometry,
                            options.chord_height,
                            max_angle,
                        );
                        // Crease flips (see `crease_flip`): a first pass
                        // free of the spread bound, refinement to restore
                        // it, then a pass that keeps it. A face the first
                        // pass leaves alone is done.
                        let geometry = &sof.geometry;
                        let chord = options.chord_height;
                        let flipped = crate::crease_flip::flip_creases(
                            &mut fm,
                            geometry,
                            chord,
                            max_angle,
                            f64::INFINITY,
                        );
                        if !flipped.is_empty() {
                            crate::angular_refine::refine_face_mesh_from(
                                &mut fm, geometry, max_angle, &flipped,
                            );
                            crate::crease_flip::flip_creases(
                                &mut fm, geometry, chord, max_angle, max_angle,
                            );
                        }
                        // The facet-tilt passes run again last: the crease
                        // flips judge creases, not lean, and can undo them.
                        crate::tilt_refine::finish_tilt(&mut fm, geometry, chord, max_angle);
                    }
                }
                return FacePassOutcome::BoundaryFirst(fm);
            }
            if trace {
                let kind = brep
                    .faces
                    .get(face_handle)
                    .and_then(|f| f.surface.as_ref())
                    .map(|sof| match &sof.geometry {
                        popkern_core::geometry::SurfaceGeometry::Plane { .. } => "Plane",
                        popkern_core::geometry::SurfaceGeometry::Cylinder { .. } => "Cylinder",
                        popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "Nurbs",
                        popkern_core::geometry::SurfaceGeometry::Gregory(_) => "Gregory",
                        _ => "Other",
                    })
                    .unwrap_or("None");
                eprintln!("  [bf] face {face_handle:?} ({kind}) fell back");
            }
            let mut mesh = crate::tessellate::tessellate_face_with_samples(
                brep,
                face_handle,
                options,
                Some(edge_samples),
            );
            if let (Ok(fm), Some(sof)) = (
                mesh.as_mut(),
                brep.faces.get(face_handle).and_then(|f| f.surface.as_ref()),
            ) {
                crate::vertex_normal::repair_singular_normals(fm, &sof.geometry);
            }
            FacePassOutcome::UniformGrid {
                seam_spine: crate::tessellate::outer_loop_is_seam_spine(brep, face_handle),
                mesh,
            }
        };
    dispatch(items, trace, per_face)
}

/// Faces incident to any edge in `edges` (both half-edge sides).
///
/// Drives the incremental encroachment rounds: a face whose boundary
/// samples did not change is a pure function of unchanged inputs, so
/// re-running it can only reproduce requests that were already
/// inserted or permanently deduped (fractions only grow, per-edge caps
/// only fill) — only faces returned here need re-checking.
pub(crate) fn faces_incident_to_edges(
    brep: &BRep,
    edges: &BTreeSet<EdgeHandle>,
) -> HashSet<FaceHandle> {
    let mut faces = HashSet::new();
    for &eh in edges {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        for h in [edge.half_edge, he.twin] {
            if let Some(hd) = brep.half_edges.get(h)
                && let Some(lp) = brep.loops.get(hd.loop_)
            {
                faces.insert(lp.face);
            }
        }
    }
    faces
}

/// Map `per_face` over `items`, parallel when the `parallel` feature
/// (default) is on and tracing is off. The parallel arm uses an
/// indexed `map` + `collect`, which rayon guarantees to be positionally
/// ordered — output is identical to the sequential arm.
pub(crate) fn dispatch<I: Send, T: Send>(
    items: Vec<I>,
    trace: bool,
    per_face: impl Fn(I) -> T + Sync,
) -> Vec<T> {
    #[cfg(feature = "parallel")]
    {
        if !trace {
            use rayon::prelude::*;
            return items.into_par_iter().map(&per_face).collect();
        }
    }
    #[cfg(not(feature = "parallel"))]
    let _ = trace;
    items.into_iter().map(per_face).collect()
}

#[cfg(test)]
mod tests {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    /// Every edge in the set contributes both of its half-edge sides.
    #[test]
    fn faces_incident_to_edges_returns_both_sides() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        for eh in brep.edges.keys() {
            let mut set = std::collections::BTreeSet::new();
            set.insert(eh);
            let faces = super::faces_incident_to_edges(&brep, &set);
            assert_eq!(
                faces.len(),
                2,
                "cube edge {eh:?} must report exactly its two incident faces"
            );
        }
    }

    /// PK-277 acceptance: identical tessellation 1-thread vs N-thread.
    /// The cylinder exercises closed rims, the periodic wall chart, and
    /// the boundary-first/uniform-grid split; meshes must match bitwise.
    #[cfg(feature = "parallel")]
    #[test]
    fn thread_count_does_not_change_mesh() {
        let brep = BRep::cylinder(0.35, 1.0, Tolerance::default()).unwrap();
        let options = crate::TessellationOptions::default();

        let single = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(|| crate::tessellate(&brep, &options).unwrap());
        let many = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .unwrap()
            .install(|| crate::tessellate(&brep, &options).unwrap());

        assert_eq!(single.triangles, many.triangles, "triangle lists differ");
        assert_eq!(
            single.vertices.len(),
            many.vertices.len(),
            "vertex counts differ"
        );
        for (i, (l, r)) in single.vertices.iter().zip(many.vertices.iter()).enumerate() {
            let bits =
                |p: &popkern_core::geometry::Point3| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
            assert_eq!(
                bits(&l.position),
                bits(&r.position),
                "vertex {i} position differs"
            );
        }
    }
}
