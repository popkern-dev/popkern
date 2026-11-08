//! The measured assignment of a low-valence PnS cap's sub-faces to its
//! patches. Split from `pns_cap_authoring.rs`; see PK-192's ledger for
//! the measured provenance.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

/// Measure which patch belongs on which sub-face: for every
/// (face, patch) pair, the worst nearest-point distance from the
/// face's rim edge curves to the patch's `v = 1` rim row. Greedy
/// bijective assignment over ascending distance; falls back to the
/// positional identity for anything left unpaired (degenerate input).
pub(crate) fn measure_sub_face_patch_assignment(
    brep: &BRep,
    sub_faces: &[FaceHandle],
    surfaces: &[SurfaceGeometry],
) -> Vec<usize> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let n = sub_faces.len();
    let mut scores: Vec<(f64, usize, usize)> = Vec::with_capacity(n * n);
    for (fi, &face) in sub_faces.iter().enumerate() {
        // Collect this face's rim sample points (edges whose curve is
        // a materialized boundary segment — i.e., every outer-loop
        // edge that predates the poke and still carries a curve).
        let mut rim_samples: Vec<popkern_core::geometry::Point3> = Vec::new();
        if let Some(ld) = brep
            .faces
            .get(face)
            .and_then(|f| brep.loops.get(f.outer_loop))
        {
            let start = ld.half_edge;
            let mut he = start;
            for _ in 0..8 {
                let Some(hd) = brep.half_edges.get(he) else {
                    break;
                };
                let next = hd.next;
                // Poke-created spoke edges carry Line placeholders to
                // the central vertex; rim edges carry the materialized
                // segment curves. Distinguish by whether the twin's
                // face is another cap sub-face-to-be — cheaper and
                // orientation-free: skip edges shorter than nothing;
                // simply sample every edge, spokes included, but only
                // rim samples separate the candidates (spokes meet at
                // the shared pole for every patch).
                if let Some(curve) = brep.edges.get(hd.edge).and_then(|e| e.curve.as_ref()) {
                    let g = &curve.geometry;
                    for k in 1..8 {
                        let t = curve.t_start + (curve.t_end - curve.t_start) * f64::from(k) / 8.0;
                        rim_samples.push((g as &dyn CurveEval).point_at(t));
                    }
                }
                he = next;
                if he == start {
                    break;
                }
            }
        }
        for (si, surface) in surfaces.iter().enumerate() {
            let mut worst = 0.0_f64;
            for p in &rim_samples {
                let mut nearest = f64::MAX;
                for j in 0..=32 {
                    let u = f64::from(j) / 32.0;
                    for v in [1.0, 0.0] {
                        nearest = nearest.min(surface.point_at(u, v).distance_to(p));
                    }
                    // Spoke isos too, so spoke samples don't penalise
                    // the true patch.
                    nearest = nearest.min(surface.point_at(0.0, u).distance_to(p));
                }
                worst = worst.max(nearest);
            }
            scores.push((worst, fi, si));
        }
    }
    scores.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut assignment: Vec<Option<usize>> = vec![None; n];
    let mut surface_used = vec![false; n];
    for (_, fi, si) in scores {
        if assignment[fi].is_none() && !surface_used[si] {
            assignment[fi] = Some(si);
            surface_used[si] = true;
        }
    }
    assignment
        .into_iter()
        .enumerate()
        .map(|(i, a)| a.unwrap_or(i))
        .collect()
}
