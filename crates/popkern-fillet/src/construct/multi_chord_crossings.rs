//! Pairwise chord-chord crossing precompute for the multi-chord
//! face split. Split from `multi_chord.rs` (file-size ratchet).

use super::multi_chord_span_geometry::{intersect_segments, intersect_trim_curves};
use super::phase_split::FaceTrim;
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Per-trim sorted crossing lists `(param_on_chord, point,
/// other_trim_index)`. Straight anchor-chord intersections first;
/// curved-support crossings fall back to the sampled trim curves
/// (PK-141), SPAN-GATED to both anchor chords: trims are authored
/// with overhang beyond their fillet's actual extent (padded offset
/// intersection), so two trims can cross far outside both strips —
/// measured: the skewed-rib junction rail's overhang crossed a rim
/// rail at x = 38.5, 4.6 beyond the rib end, minting a phantom
/// intermediate vertex whose sub-chord straddled two loops (PK-185).
/// A real crossing lies within both anchor chords' axial spans
/// (curved-support crossings sit near, not on, their chords — the
/// gate is on the PROJECTION parameter, not the distance).
pub(super) fn precompute_chord_crossings(
    brep: &BRep,
    trims: &[FaceTrim],
    tolerance: &Tolerance,
) -> Vec<Vec<(f64, Point3, usize)>> {
    let n = trims.len();
    let segments: Vec<(Point3, Point3)> = trims
        .iter()
        .map(|t| {
            (
                brep.vertices[t.anchors.0].position,
                brep.vertices[t.anchors.1].position,
            )
        })
        .collect();
    let mut intersections: Vec<Vec<(f64, Point3, usize)>> = vec![Vec::new(); n];
    let cross_gate = tolerance.linear().max(1e-9) * 1.0e4;
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = segments[i];
            let (c, d) = segments[j];
            if let Some((s, t, pt)) = intersect_segments(a, b, c, d) {
                intersections[i].push((s, pt, j));
                intersections[j].push((t, pt, i));
            } else if let Some((s, t, pt)) =
                intersect_trim_curves(&trims[i].trim_3d, &trims[j].trim_3d, cross_gate)
            {
                let within_span = |p0: Point3, p1: Point3| -> bool {
                    let dseg = p1 - p0;
                    let len_sq = dseg.dot(&dseg);
                    if len_sq < 1e-18 {
                        return false;
                    }
                    let f = (pt - p0).dot(&dseg) / len_sq;
                    (-0.02..=1.02).contains(&f)
                };
                if within_span(a, b) && within_span(c, d) {
                    intersections[i].push((s, pt, j));
                    intersections[j].push((t, pt, i));
                }
            }
        }
    }
    for list in intersections.iter_mut() {
        list.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    intersections
}
