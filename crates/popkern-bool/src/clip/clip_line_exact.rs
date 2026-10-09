//! Exact line clipping against a planar face's TRUE boundary curves
//! (PK-751, PK-753).
//!
//! The planar split fallback extends a trim line across its face and
//! clips it to the face boundary, so each fragment gets a cut that ends
//! on its own edges. Its Cyrus–Beck clipper works on the outer loop's
//! VERTEX polygon, which on a face bounded by arcs is a chord polygon:
//! the cut's ends landed on chords, up to a sagitta inside the true
//! boundary, at points that are no vertex and on no edge, and the face
//! split refused ("trim curve start point does not match any vertex";
//! measured 5.8e-4, 0.51 and 3.1 on the six-axis arm's castings, where
//! a web's plane crosses a ring's flat face out to its rim). Cyrus–Beck
//! is also convex-only.
//!
//! This clip intersects the line with every boundary edge curve of the
//! face (outer and inner loops): sign changes of the in-plane offset
//! over a sampled edge, each root refined by bisection on the edge's own
//! parameter. Crossings sorted along the line alternate out/in; the
//! span holding the trim's midpoint is the clip, and its ends lie on
//! the boundary curves to bisection precision.
//!
//! Offsets are classified half-open (zero counts as negative) per edge
//! in the edge's own direction, so a line through a shared vertex is
//! counted once, as in the classic ray-crossing rule.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

/// Samples per edge for the sign-change scan. A crossing pair closer
/// together than an edge's 1/64 is a near-tangency the scan may miss;
/// missing both crossings of a pair keeps the parity.
const EDGE_SAMPLES: usize = 64;

/// Bisection steps per root: 2^-60 of the bracketing sample span.
const BISECTIONS: usize = 60;

/// Clip the line `line` (degree 1, parameter `[0, 1]`) to the span of
/// the planar face `face` that holds the line parameter `s_mid`.
/// `None` when the face is not planar, `s_mid` lies outside the face,
/// or no boundary crossing brackets it.
pub(crate) fn clip_line_to_planar_face_exact(
    brep: &BRep,
    face: FaceHandle,
    line: &NurbsCurve,
    s_mid: f64,
) -> Option<(f64, f64)> {
    let face_data = brep.faces.get(face)?;
    let SurfaceGeometry::Plane { normal, .. } = &face_data.surface.as_ref()?.geometry else {
        return None;
    };
    let n = normal.normalized()?;
    let (t0, t1) = line.parameter_range();
    let a = line.point_at(t0);
    let d = line.point_at(t1) - a;
    let d2 = d.dot(&d);
    if d2 <= 0.0 {
        return None;
    }
    // In-plane offset of q from the line (scaled by |d|).
    let offset = |q: &Point3| (*q - a).cross(&d).dot(&n);
    let along = |q: &Point3| (*q - a).dot(&d) / d2;

    let mut crossings: Vec<f64> = Vec::new();
    let loops = std::iter::once(face_data.outer_loop).chain(face_data.inner_loops.iter().copied());
    for lh in loops {
        let start = brep.loops.get(lh)?.half_edge;
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let h = brep.half_edges.get(he)?;
            let curve = brep.edges.get(h.edge)?.curve.as_ref()?;
            // Walk the edge in the half-edge's direction (the shared
            // vertex convention needs a consistent traversal).
            let origin = brep.vertices.get(h.origin)?.position;
            let forward = curve.geometry.point_at(curve.t_start).distance_to(&origin)
                <= curve.geometry.point_at(curve.t_end).distance_to(&origin);
            let (ta, tb) = if forward {
                (curve.t_start, curve.t_end)
            } else {
                (curve.t_end, curve.t_start)
            };
            let at = |k: usize| ta + (tb - ta) * k as f64 / EDGE_SAMPLES as f64;
            let positive = |t: f64| offset(&curve.geometry.point_at(t)) > 0.0;
            let mut prev = positive(at(0));
            for k in 1..=EDGE_SAMPLES {
                let now = positive(at(k));
                if now != prev {
                    let (mut lo, mut hi) = (at(k - 1), at(k));
                    for _ in 0..BISECTIONS {
                        let mid = 0.5 * (lo + hi);
                        if positive(mid) == prev {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    crossings.push(along(&curve.geometry.point_at(0.5 * (lo + hi))));
                }
                prev = now;
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    crossings.sort_by(f64::total_cmp);
    // Out/in alternation: the midpoint is inside after an odd count.
    let before = crossings.iter().filter(|&&s| s < s_mid).count();
    if before % 2 == 0 {
        return None;
    }
    let enter = crossings[before - 1];
    let exit = *crossings.get(before)?;
    Some((enter, exit))
}

#[cfg(test)]
#[path = "clip_line_exact_tests.rs"]
mod tests;
