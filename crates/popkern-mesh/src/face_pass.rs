//! Per-face tessellation pass dispatch (PK-277).
//!
//! Both per-face passes of `tessellate_brep` — encroachment
//! split-request collection and final face meshing — are pure
//! functions of `&BRep` (plain `SlotMap` storage, no interior
//! mutability), the shared read-only canonical edge samples, and the
//! options, so faces can run on separate threads. Per-face results are
//! collected positionally (rayon's indexed `collect`) and flattened in
//! caller face order, so output is identical at any thread count.
//! `PK_BF_TRACE` forces the sequential path so interleaved thread
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
    /// The face fell back to the legacy pcurve-grid path.
    Legacy {
        /// True when the face's outer loop is a pure seam spine
        /// (every edge traversed twice) — such faces keep their exact
        /// chart sampling and are excluded from projection sewing.
        seam_spine: bool,
        /// The legacy tessellation result; the caller propagates the
        /// first error in face order.
        mesh: Result<FaceMesh, MeshError>,
    },
}

/// Collect encroachment split requests for `faces`, in face order.
///
/// Each face's request list is computed independently; the flattened
/// sequence is identical to running the faces sequentially, so the
/// caller's insert-and-dedup pass sees the same requests in the same
/// order at any thread count.
pub(crate) fn collect_split_requests(
    brep: &BRep,
    faces: &[FaceHandle],
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
    options: &TessellationOptions,
    trace: bool,
) -> Vec<SplitRequest> {
    let per_face = |face_handle: FaceHandle| -> Vec<SplitRequest> {
        let mut requests = Vec::new();
        crate::boundary_first::tessellate_face_boundary_first(
            brep,
            face_handle,
            edge_samples,
            options,
            &mut requests,
        );
        requests
    };
    let nested = dispatch(faces, trace, per_face);
    nested.into_iter().flatten().collect()
}

/// Run the final tessellation pass for `faces`, in face order.
///
/// Boundary-first is attempted per face; fallbacks run the legacy path
/// with canonical samples. Legacy errors are carried in the outcome so
/// the caller can propagate the first one in face order — positionally
/// identical to the retired sequential loop's early `?` return.
pub(crate) fn tessellate_face_pass(
    brep: &BRep,
    faces: &[FaceHandle],
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
    options: &TessellationOptions,
    trace: bool,
) -> Vec<FacePassOutcome> {
    let per_face = |face_handle: FaceHandle| -> FacePassOutcome {
        let mut sink = Vec::new();
        if let Some(fm) = crate::boundary_first::tessellate_face_boundary_first(
            brep,
            face_handle,
            edge_samples,
            options,
            &mut sink,
        ) {
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
        FacePassOutcome::Legacy {
            seam_spine: crate::tessellate::outer_loop_is_seam_spine(brep, face_handle),
            mesh: crate::tessellate::tessellate_face_with_samples(
                brep,
                face_handle,
                options,
                Some(edge_samples),
            ),
        }
    };
    dispatch(faces, trace, per_face)
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

/// Map `per_face` over `faces`, parallel when the `parallel` feature
/// (default) is on and tracing is off. The parallel arm uses an
/// indexed `map` + `collect`, which rayon guarantees to be positionally
/// ordered — output is identical to the sequential arm.
fn dispatch<T: Send>(
    faces: &[FaceHandle],
    trace: bool,
    per_face: impl Fn(FaceHandle) -> T + Sync,
) -> Vec<T> {
    #[cfg(feature = "parallel")]
    {
        if !trace {
            use rayon::prelude::*;
            return faces.par_iter().map(|&fh| per_face(fh)).collect();
        }
    }
    #[cfg(not(feature = "parallel"))]
    let _ = trace;
    faces.iter().map(|&fh| per_face(fh)).collect()
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
    /// the boundary-first/legacy split; meshes must match bitwise.
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
