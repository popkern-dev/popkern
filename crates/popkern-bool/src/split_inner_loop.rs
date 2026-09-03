//! Inner-loop detection helpers for the split stage.
//!
//! When a closed intersection curve has 0 boundary crossings and lies
//! entirely inside a face, it should become an inner loop (hole) rather
//! than being dropped.

use popkern_core::{
    geometry::Point3,
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::clip::face_uv_boundary;

/// Test whether a 2D point (using x,y of a `Point3`) is inside a polygon
/// via the even-odd ray-casting rule.
///
/// Casts a horizontal ray from `(px, py)` in the +x direction and counts
/// crossings with the polygon edges.
pub(crate) fn ray_cast_contains(polygon: &[Point3], px: f64, py: f64) -> bool {
    let n = polygon.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let yi = polygon[i].y;
        let yj = polygon[j].y;
        if (yi > py) != (yj > py) {
            let xi = polygon[i].x;
            let xj = polygon[j].x;
            let x_intersect = xi + (py - yi) / (yj - yi) * (xj - xi);
            if px < x_intersect {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// Check whether a UV-space point lies inside a face's outer boundary.
///
/// Samples the face boundary pcurves to build a polygon, then applies
/// even-odd ray casting for containment.
pub(crate) fn point_inside_face_uv(
    brep: &BRep,
    face: FaceHandle,
    uv_point: Point3,
    tolerance: &Tolerance,
) -> bool {
    let boundary = match face_uv_boundary(brep, face) {
        Some(b) => b,
        None => return false,
    };
    // A face whose EVERY boundary edge is used twice by its own loops
    // (a revolved circle's torus: its tube seam and its equator seam,
    // nothing else) has no genuine boundary — it covers its whole
    // chart, and every on-surface query is inside. Neither periodic
    // route below can say so: the chart cell declines a ring that
    // winds both axes, and the u-parity count never folds v, so a
    // box subtraction found both of its intersection circles
    // "outside" and served the torus untouched (PK-418).
    if crate::clip::face_covers_full_chart(brep, face) {
        return true;
    }
    // SEAM-BRANCH chaining: on a periodic chart a seam edge carries
    // one pcurve shared by both of its uses, so the raw contour
    // jumps a full period mid-loop and even-odd containment reads
    // every interior point as OUTSIDE (measured: the recess
    // cylinder's mouth ellipse — a closed record wholly interior to
    // the wall chart — probed inside=false, the record was dropped,
    // and the seat subtract degraded to the no-records containment
    // path, leaving an untrimmed chimney in the block). Rebranch a
    // pcurve ONLY when its start sits an EXACT whole-period multiple
    // from the previous end: a seam-branch offset is exact by
    // construction, while the arbitrary-u jumps at a cone apex are
    // not — blind chaining folded the cone wall's contour and its
    // through-circles were dropped instead (cone/sphere inner-loop
    // suites measured the distinction).
    let periodicity = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.periodicity())
        .unwrap_or((None, None));
    let mut polygon = Vec::new();
    // Per-pcurve sample RUNS (segments only within one pcurve — the
    // even-odd polygon's implicit closure between unrelated boundary
    // pieces fabricated phantom edges on rotational charts).
    let mut runs: Vec<Vec<Point3>> = Vec::new();
    let sample_contour = |contour: &[(
        popkern_core::topology::HalfEdgeHandle,
        popkern_core::geometry::pcurve::Pcurve,
    )],
                          samples_per_edge: usize,
                          polygon: &mut Vec<Point3>,
                          runs: &mut Vec<Vec<Point3>>| {
        for (_he, pcurve) in contour {
            let (t0, t1) = pcurve.parameter_range();
            let range = t1 - t0;
            if range.abs() < tolerance.linear() {
                continue;
            }
            let mut run = Vec::with_capacity(samples_per_edge + 1);
            for i in 0..=samples_per_edge {
                let t = t0 + range * (i as f64) / (samples_per_edge as f64);
                let q = pcurve.point_at(t);
                run.push(Point3::new(q.x, q.y, 0.0));
                // The planar even-odd polygon keeps the historical
                // endpoint-exclusive sampling — the tuned verdicts of
                // marginal planar cases depend on it.
                if i < samples_per_edge {
                    polygon.push(Point3::new(q.x, q.y, 0.0));
                }
            }
            runs.push(run);
        }
    };
    sample_contour(&boundary, 8, &mut polygon, &mut runs);
    if periodicity.0.is_none() && periodicity.1.is_none() {
        return ray_cast_contains(&polygon, uv_point.x, uv_point.y);
    }
    // On a U-PERIODIC chart, count parity along the AXIAL direction
    // using only REAL boundary segments: intersect each pcurve's own
    // sample runs with the vertical line u = u_q and count crossings
    // with v above the query. The even-odd polygon's +U ray was
    // hostage to the seam's period branch, and its implicit closure
    // segments between unrelated boundary pieces (seam uses, apex
    // degenerates) fabricated phantom edges — a full-period wall's
    // interior read OUTSIDE and its through-records were dropped
    // (measured on the recess mouth ellipse and the cone through-
    // circles). Seam edges are constant-U and can never cross the
    // axial line; rims cross it exactly once per period image.
    let u_probe: Vec<f64> = match periodicity.0 {
        Some(p) => vec![uv_point.x, uv_point.x + p, uv_point.x - p],
        None => vec![uv_point.x],
    };
    let mut crossings_above = 0usize;
    let mut has_u_extent = false;
    let u_span_floor = periodicity.0.map(|p| p * 1e-9).unwrap_or(1e-12);
    for run in &runs {
        for w in run.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            if (a.x - b.x).abs() > u_span_floor {
                has_u_extent = true;
            }
            for &uq in &u_probe {
                if (a.x > uq) != (b.x > uq) {
                    let f = (uq - a.x) / (b.x - a.x);
                    let v = a.y + f * (b.y - a.y);
                    if v > uv_point.y {
                        crossings_above += 1;
                    }
                }
            }
        }
    }
    if !has_u_extent {
        // A boundary made only of seam-vertical and degenerate
        // pieces (a full primitive sphere: seam + poles) excludes
        // nothing — the face covers the whole chart, and every
        // on-surface query is inside.
        return true;
    }
    crossings_above % 2 == 1
}

/// Hole-aware REGION containment: [`point_inside_face_uv`] minus the
/// face's holes (even-odd XOR per inner loop). The tuned crossing and
/// split machinery deliberately keeps the OUTER-ONLY test — records
/// running THROUGH holes are partitioned at split time by the
/// PK-185/PK-186 inner-loop band machinery, and clipping them at the
/// clip stage stranded three tuned countersink/socket fixtures
/// (measured 2026-09-03). This variant serves the two callers that
/// mean the true face REGION: the zero-crossing clip verdict (a
/// record wholly inside a hole is off the face — the vent-slot class)
/// and the island-inside check (an island accepted inside a hole
/// imprinted an invalid imbricated loop — the lid∪lip ring class).
/// Hole contours sample at 32 per pcurve: a full-circle hole at 8
/// samples is an octagon whose inradius sits 1.4 % under the rim and
/// near-rim points misclassify.
/// [`point_inside_face_region_uv`] scoped to PRE-EXISTING holes: a
/// hole whose every edge is a freshly-minted intersection edge of the
/// running boolean is skipped. An island landing inside such a fresh
/// hole is a legitimately NESTED imprint (the annular flush contact
/// of a ring lip: outer outline first — a hole plus its enclosed
/// piece face — then the inner outline inside it), while an island
/// inside a hole the face brought INTO the boolean is off the face
/// region (the vent-slot class; accepting it imprinted an invalidly
/// imbricated loop).
pub(crate) fn point_inside_face_region_excluding_fresh_holes(
    brep: &BRep,
    face: FaceHandle,
    uv_point: Point3,
    tolerance: &Tolerance,
    fresh_edges: &std::collections::HashSet<popkern_core::topology::EdgeHandle>,
) -> bool {
    if !point_inside_face_uv(brep, face, uv_point, tolerance) {
        return false;
    }
    let mut inside = true;
    for contour in crate::clip::face_uv_inner_boundaries(brep, face) {
        let fresh = !contour.is_empty()
            && contour.iter().all(|(he, _)| {
                brep.half_edges
                    .get(*he)
                    .is_some_and(|h| fresh_edges.contains(&h.edge))
            });
        if fresh {
            continue;
        }
        if hole_contour_contains(&contour, uv_point, tolerance) {
            inside = !inside;
        }
    }
    inside
}

/// Even-odd test of one hole contour, sampled at 32 per pcurve (an
/// 8-sample circle is an octagon whose inradius sits 1.4 % under the
/// rim; near-rim points misclassify).
fn hole_contour_contains(
    contour: &[(
        popkern_core::topology::HalfEdgeHandle,
        popkern_core::geometry::pcurve::Pcurve,
    )],
    uv_point: Point3,
    tolerance: &Tolerance,
) -> bool {
    let mut hole_polygon = Vec::new();
    let samples_per_edge = 32;
    for (_he, pcurve) in contour {
        let (t0, t1) = pcurve.parameter_range();
        let range = t1 - t0;
        if range.abs() < tolerance.linear() {
            continue;
        }
        for i in 0..samples_per_edge {
            let t = t0 + range * (i as f64) / (samples_per_edge as f64);
            let q = pcurve.point_at(t);
            hole_polygon.push(Point3::new(q.x, q.y, 0.0));
        }
    }
    ray_cast_contains(&hole_polygon, uv_point.x, uv_point.y)
}

pub(crate) fn point_inside_face_region_uv(
    brep: &BRep,
    face: FaceHandle,
    uv_point: Point3,
    tolerance: &Tolerance,
) -> bool {
    if !point_inside_face_uv(brep, face, uv_point, tolerance) {
        return false;
    }
    let mut inside = true;
    for contour in crate::clip::face_uv_inner_boundaries(brep, face) {
        if hole_contour_contains(&contour, uv_point, tolerance) {
            inside = !inside;
        }
    }
    inside
}

#[cfg(test)]
#[path = "split_inner_loop_tests.rs"]
mod tests;
