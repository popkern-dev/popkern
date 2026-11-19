//! The atlas's lattice rung (chart-service.md, step 7: the grid fallback
//! on the atlas). A face reaching a chart singularity whose structured
//! atlas meshing declined on the final edge samples meshes again on the
//! same atlas with no refinement rounds: a crossing of a cap row or of the
//! band's slit meridian, which the structured rung requests as a canonical
//! sample, is inserted into the ring where it lies on its edge, and the
//! band meshes as a quadtree lattice in the face chart clipped by its
//! rings. The pieces still lie in charts with no singular row, so this
//! rung has no pole handling of its own: one authority for poles in every
//! mesh path.

use super::*;
#[path = "boundary_first_atlas_lattice_band.rs"]
mod band;
pub(super) use band::mesh_band;

/// Insertion rounds before the lattice rung declines: each inserts every
/// crossing the last round asked for, so a face whose rings cross a row
/// and the slit needs two.
const ROUNDS: usize = 6;

/// Mesh `face_handle` on the atlas's lattice rung: `None` when the face
/// reaches no cap (the atlas does not own it), `Some(None)` when the rung
/// declined (its rings could not be built or assembled into pieces).
pub(in crate::boundary_first) fn tessellate_on_atlas_lattice(
    brep: &BRep,
    face_handle: FaceHandle,
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
    options: &TessellationOptions,
) -> Option<Option<FaceMesh>> {
    let face = brep.faces.get(face_handle)?;
    let surface_on_face = face.surface.as_ref()?;
    let (surface, orientation) = (&surface_on_face.geometry, surface_on_face.orientation);
    let descriptor = ChartDescriptor::of(surface, &brep.tolerance);
    if !descriptor
        .singular_rows
        .iter()
        .any(|r| r.kind == SingularKind::Chart)
    {
        return None;
    }
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());
    // A loop with no edges (a vertex loop, a STEP face bounded by one
    // alone) bounds nothing: no ring. A face with only such loops is the
    // whole surface, its pieces bounded by the cap rows alone.
    let built: Option<Vec<BoundaryRing>> = loops
        .into_iter()
        .filter(|&lh| {
            brep.loops
                .get(lh)
                .is_some_and(|l| brep.half_edges.contains_key(l.half_edge))
        })
        .map(|lh| super::super::build_ring(brep, lh, edge_samples))
        .collect();
    let Some(mut built) = built else {
        why(face_handle, "a boundary ring could not be built");
        return Some(None);
    };
    for _ in 0..ROUNDS {
        let mut requests = Vec::new();
        let out = super::tessellate_on_atlas(
            brep,
            face_handle,
            (surface, orientation),
            &built,
            (&descriptor, options, PieceMeshing::Lattice),
            &mut requests,
        )?;
        if let Some(mesh) = out {
            return Some(Some(mesh));
        }
        // Long-segment splits the cap asks for are refinement, not
        // crossings: only a decline's requests are inserted.
        let inserted = insert_samples(brep, (surface, &descriptor), &mut built, &requests);
        if inserted == 0 {
            why(face_handle, "pieces declined with no crossing to insert");
            return Some(None);
        }
    }
    why(face_handle, "crossings still requested after every round");
    Some(None)
}

fn why(face: FaceHandle, reason: &str) {
    if popkern_core::trace::trace_enabled("atlas") {
        eprintln!("    [atlas] face {face:?}: lattice rung declined: {reason}");
    }
}

/// Insert each requested `(edge, fraction)` into the ring segments whose
/// span on that edge holds it strictly: the point on the edge's curve, its
/// face parameters the surface's at that point on the segment's chart
/// branch (the periodic offset nearest the segment's own parameters).
/// Returns how many samples went in.
fn insert_samples(
    brep: &BRep,
    (surface, descriptor): (&SurfaceGeometry, &ChartDescriptor),
    rings: &mut [BoundaryRing],
    requests: &[SplitRequest],
) -> usize {
    let mut inserted = 0;
    for &(edge, f) in requests {
        if edge == EdgeHandle::default() {
            continue;
        }
        let Some(curve) = brep.edges.get(edge).and_then(|e| e.curve.as_ref()) else {
            continue;
        };
        for ring in rings.iter_mut() {
            let n = ring.uv.len();
            let Some(i) = (0..n).find(|&i| {
                let (e, fa, fb) = ring.sources[i];
                e == edge && (f - fa) * (f - fb) < 0.0
            }) else {
                continue;
            };
            let (_, fa, fb) = ring.sources[i];
            let j = (i + 1) % n;
            let p = curve
                .geometry
                .point_at(curve.t_start + (curve.t_end - curve.t_start) * f);
            let s = (f - fa) / (fb - fa);
            let (a, b) = (ring.uv[i], ring.uv[j]);
            let guess = (a.0 + (b.0 - a.0) * s, a.1 + (b.1 - a.1) * s);
            let uv = surface
                .closest_point(&p)
                .map(|q| on_branch_of(q, guess, descriptor))
                .unwrap_or(guess);
            ring.uv.insert(i + 1, uv);
            ring.pos.insert(i + 1, p);
            ring.sources[i] = (edge, fa, f);
            ring.sources.insert(i + 1, (edge, f, fb));
            inserted += 1;
        }
    }
    inserted
}

/// `q` moved by whole periods along each periodic axis to the branch
/// nearest `guess`.
pub(super) fn on_branch_of(
    q: (f64, f64),
    guess: (f64, f64),
    descriptor: &ChartDescriptor,
) -> (f64, f64) {
    use popkern_core::geometry::chart::ops::nearest_offset;
    let shift = |x: f64, g: f64, axis: Axis| match descriptor.period(axis) {
        Some(p) => g + nearest_offset(x, g, p.period),
        None => x,
    };
    (shift(q.0, guess.0, Axis::U), shift(q.1, guess.1, Axis::V))
}
