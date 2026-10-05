//! Adaptive trim polyline construction from face boundary loops.
//!
//! Builds UV-space polylines from face loops using pcurves (direct evaluation)
//! with 3D-to-UV projection fallback. These polylines are used by later tasks
//! for triangle classification and clipping.

// Items are pub(crate) for use by the trim pipeline modules.

use crate::trim_sampling::sample_loop_with_samples;
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_math::PortableMath;

/// A closed polyline in UV space representing a face boundary loop.
pub(crate) struct TrimLoop {
    /// UV points along the boundary.
    pub points: Vec<(f64, f64)>,
    /// 3D positions corresponding to each UV point, used for validating
    /// clipping vertex positions when UV projection is imprecise.
    pub points_3d: Vec<Point3>,
    /// Whether this is the outer boundary (`true`) or a hole (`false`).
    pub is_outer: bool,
    /// Whether this loop was built entirely from pcurves (reliable) or
    /// used 3D→UV projection fallback (may be unreliable for non-planar
    /// surfaces).
    #[allow(dead_code)]
    pub from_pcurves: bool,
}

/// Build trim loops for all boundary loops of a face.
///
/// Walks the outer loop and all inner loops, sampling each into a UV polyline.
/// The outer loop is tagged with `is_outer: true`; inner loops with `false`.
/// Loops with fewer than 3 points are skipped.
/// Build trim loops; canonical per-edge stations when provided
/// (PK-299 — see `trim_canonical`).
pub(crate) fn build_trim_loops_with_samples(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    chord_height: f64,
    edge_samples: Option<
        &std::collections::HashMap<
            popkern_core::topology::EdgeHandle,
            crate::boundary_first::EdgeSamples,
        >,
    >,
) -> Vec<TrimLoop> {
    let Some(face) = brep.faces.get(face_handle) else {
        return Vec::new();
    };

    let (period_u, period_v) = surface.periodicity();
    let mut result = Vec::new();

    // Outer loop.
    let (mut outer_points, mut outer_points_3d, outer_from_pcurves) =
        sample_loop_with_samples(brep, face.outer_loop, surface, chord_height, edge_samples);
    if outer_points.len() >= 3 {
        // Unwrap seam-crossing polygons on periodic surfaces.
        if crosses_seam(&outer_points, period_u, period_v) {
            unwrap_periodic(&mut outer_points, period_u, period_v);
        }
        // Canonical winding: outer loops CCW in UV. Half-edge order
        // around a face is CCW in 3D but its UV image winding depends
        // on the surface parameterization's handedness (measured:
        // structurally identical fillet faces come out with opposite
        // signs), and the inside tests assume CCW outers.
        dedup_loop(&mut outer_points, &mut outer_points_3d);
        normalize_winding(&mut outer_points, &mut outer_points_3d, true);
        result.push(TrimLoop {
            points: outer_points,
            points_3d: outer_points_3d,
            is_outer: true,
            from_pcurves: outer_from_pcurves,
        });
    }

    // Inner loops (holes): canonical winding CW.
    for &inner_loop in &face.inner_loops {
        let (mut inner_points, mut inner_points_3d, inner_from_pcurves) =
            sample_loop_with_samples(brep, inner_loop, surface, chord_height, edge_samples);
        if inner_points.len() >= 3 {
            if crosses_seam(&inner_points, period_u, period_v) {
                unwrap_periodic(&mut inner_points, period_u, period_v);
            }
            dedup_loop(&mut inner_points, &mut inner_points_3d);
            normalize_winding(&mut inner_points, &mut inner_points_3d, false);
            result.push(TrimLoop {
                points: inner_points,
                points_3d: inner_points_3d,
                is_outer: false,
                from_pcurves: inner_from_pcurves,
            });
        }
    }

    result
}

/// Remove consecutive duplicate points (each half-edge contributes
/// its start AND end, doubling shared corners — the zero-length
/// polygon edges break ray-cast parity along their v-rows) and the
/// closing duplicate.
fn dedup_loop(points: &mut Vec<(f64, f64)>, points_3d: &mut Vec<Point3>) {
    let mut i = 1;
    while i < points.len() {
        let (px, py) = points[i - 1];
        let (qx, qy) = points[i];
        if (px - qx).abs() < 1.0e-12 && (py - qy).abs() < 1.0e-12 {
            points.remove(i);
            if i < points_3d.len() {
                points_3d.remove(i);
            }
        } else {
            i += 1;
        }
    }
    if points.len() >= 2 {
        let first = points[0];
        let last = points[points.len() - 1];
        if (first.0 - last.0).abs() < 1.0e-12 && (first.1 - last.1).abs() < 1.0e-12 {
            points.pop();
            points_3d.pop();
        }
    }
}

/// Reverse a loop in place when its signed area disagrees with the
/// requested orientation (`ccw = true` → positive area).
fn normalize_winding(points: &mut [(f64, f64)], points_3d: &mut [Point3], ccw: bool) {
    let mut area2 = 0.0;
    for i in 0..points.len() {
        let (x1, y1) = points[i];
        let (x2, y2) = points[(i + 1) % points.len()];
        area2 += x1 * y2 - x2 * y1;
    }
    if (area2 > 0.0) != ccw {
        points.reverse();
        points_3d.reverse();
    }
}

/// Sample a single loop into UV-space points with corresponding 3D positions.
///
/// For each half-edge in the loop:
/// - If the half-edge has a pcurve, samples it directly in UV space via
///   [`adaptive_sample_pcurve`], then evaluates the surface to get 3D positions.
/// - If no pcurve but the edge has a 3D curve, projects 32 uniform curve
///   samples to UV via [`project_point_to_uv`], respecting half-edge direction,
///   and stores the original 3D curve points.
/// - If neither, projects the origin vertex position to UV.
///
/// Returns `(uv_points, points_3d, all_from_pcurves)`.
/// Detect if a UV polyline crosses a periodic seam.
///
/// A seam crossing is detected when consecutive points in the polyline have a
/// gap larger than 40% of the period in either periodic direction. This
/// threshold reliably catches genuine wraps while ignoring small numerical
/// jitter.
pub(crate) fn crosses_seam(
    points: &[(f64, f64)],
    period_u: Option<f64>,
    period_v: Option<f64>,
) -> bool {
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        if let Some(pu) = period_u
            && (points[j].0 - points[i].0).abs() > pu * 0.4
        {
            return true;
        }
        if let Some(pv) = period_v
            && (points[j].1 - points[i].1).abs() > pv * 0.4
        {
            return true;
        }
    }
    false
}

/// Unwrap a UV polyline so it doesn't cross a periodic seam.
///
/// Shifts points so the polygon is contiguous in parameter space. For each
/// periodic dimension, computes the circular mean of all point coordinates,
/// then shifts each point by integer multiples of the period so it falls
/// within `period/2` of the mean.
pub(crate) fn unwrap_periodic(
    points: &mut [(f64, f64)],
    period_u: Option<f64>,
    period_v: Option<f64>,
) {
    if let Some(pu) = period_u {
        let center_u = circular_mean(points.iter().map(|p| p.0), pu);
        for p in points.iter_mut() {
            while p.0 - center_u > pu / 2.0 {
                p.0 -= pu;
            }
            while p.0 - center_u < -pu / 2.0 {
                p.0 += pu;
            }
        }
    }
    if let Some(pv) = period_v {
        let center_v = circular_mean(points.iter().map(|p| p.1), pv);
        for p in points.iter_mut() {
            while p.1 - center_v > pv / 2.0 {
                p.1 -= pv;
            }
            while p.1 - center_v < -pv / 2.0 {
                p.1 += pv;
            }
        }
    }
}

/// Circular mean of angles with a given period.
///
/// Computes the mean direction by converting each value to a unit vector on
/// the circle (with the given period mapped to 2pi), averaging the sine and
/// cosine components, and converting back.
fn circular_mean(angles: impl Iterator<Item = f64>, period: f64) -> f64 {
    let (sin_sum, cos_sum, count) = angles.fold((0.0, 0.0, 0usize), |(s, c, n), a| {
        let theta = a * std::f64::consts::TAU / period;
        (s + theta.psin(), c + theta.pcos(), n + 1)
    });
    if count == 0 {
        return 0.0;
    }
    let mean_theta = (sin_sum / count as f64).patan2(cos_sum / count as f64);
    mean_theta * period / std::f64::consts::TAU
}

/// Returns the shoelace area of the outer loop for degenerate-boundary detection.
///
/// If no outer loop is found, returns 0.0.
pub(crate) fn trim_loop_area(loops: &[TrimLoop]) -> f64 {
    for tl in loops {
        if tl.is_outer {
            return shoelace_area(&tl.points);
        }
    }
    0.0
}

/// Compute the area of a 2D polygon using the shoelace formula.
///
/// Returns the absolute value of the signed area.
pub(crate) fn shoelace_area(polygon: &[(f64, f64)]) -> f64 {
    if polygon.len() < 3 {
        return 0.0;
    }
    let mut area = 0.0_f64;
    let n = polygon.len();
    for i in 0..n {
        let (x1, y1) = polygon[i];
        let (x2, y2) = polygon[(i + 1) % n];
        area += x1 * y2 - x2 * y1;
    }
    (area / 2.0).abs()
}

/// Classification of a triangle relative to trim loops.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TriClass {
    /// Triangle lies entirely inside the trim region.
    Inside,
    /// Triangle lies entirely outside the trim region.
    Outside,
    /// Triangle straddles a trim boundary.
    Boundary,
}

/// Classify a triangle against a set of trim loops.
///
/// Computes the centroid of the triangle. If the centroid is not inside
/// the trim region, returns [`TriClass::Outside`]. If it is inside, checks
/// whether any triangle edge intersects any trim loop segment; if so,
/// returns [`TriClass::Boundary`]. Otherwise returns [`TriClass::Inside`].
pub(crate) fn classify_triangle(tri_uv: &[(f64, f64); 3], trim_loops: &[TrimLoop]) -> TriClass {
    // Boundary detection FIRST: a triangle whose edges cross the trim
    // (or which contains a trim vertex) straddles the boundary even
    // when its centroid lies outside. The previous centroid-first
    // order discarded every straddling cell with under half its area
    // inside — the systematic sawtooth along all trimmed boundaries.
    let tri_edges = [(0, 1), (1, 2), (2, 0)];
    for loop_ in trim_loops {
        let n = loop_.points.len();
        if n < 2 {
            continue;
        }
        for i in 0..n {
            let j = (i + 1) % n;
            let seg_a = loop_.points[i];
            let seg_b = loop_.points[j];
            for &(ei, ej) in &tri_edges {
                if segments_intersect(tri_uv[ei], tri_uv[ej], seg_a, seg_b) {
                    return TriClass::Boundary;
                }
            }
        }
        for &p in &loop_.points {
            if point_in_triangle_uv(p, tri_uv) {
                return TriClass::Boundary;
            }
        }
    }

    let cx = (tri_uv[0].0 + tri_uv[1].0 + tri_uv[2].0) / 3.0;
    let cy = (tri_uv[0].1 + tri_uv[1].1 + tri_uv[2].1) / 3.0;
    if point_inside_trim_loops(cx, cy, trim_loops) {
        TriClass::Inside
    } else {
        TriClass::Outside
    }
}

/// Barycentric point-in-triangle test in UV space.
pub(crate) fn point_in_triangle_uv(p: (f64, f64), tri: &[(f64, f64); 3]) -> bool {
    let sign = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
        (a.0 - c.0) * (b.1 - c.1) - (b.0 - c.0) * (a.1 - c.1)
    };
    let d1 = sign(p, tri[0], tri[1]);
    let d2 = sign(p, tri[1], tri[2]);
    let d3 = sign(p, tri[2], tri[0]);
    let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
    let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
    !(has_neg && has_pos)
}

/// Test whether a point lies inside the trim region defined by trim loops.
///
/// For outer loops a nonzero winding number means the point is inside.
/// For inner loops (holes) a nonzero winding number means the point is
/// inside the hole, so it is outside the trim region. Outer loops are
/// processed first.
pub(crate) fn point_inside_trim_loops(px: f64, py: f64, loops: &[TrimLoop]) -> bool {
    let mut inside = false;

    // Process outer loops first, then inner loops.
    for tl in loops.iter().filter(|l| l.is_outer) {
        if winding_number(px, py, &tl.points) != 0 {
            inside = true;
        }
    }

    if !inside {
        return false;
    }

    for tl in loops.iter().filter(|l| !l.is_outer) {
        if winding_number(px, py, &tl.points) != 0 {
            return false;
        }
    }

    true
}

/// Compute the winding number of a point relative to a closed polygon.
///
/// Uses the standard crossing-number algorithm: for each edge of the
/// polygon, determines whether the horizontal ray from `(px, py)` to
/// `+inf` crosses that edge, and accumulates +1 or -1 depending on the
/// crossing direction.
pub(crate) fn winding_number(px: f64, py: f64, polygon: &[(f64, f64)]) -> i32 {
    let n = polygon.len();
    if n < 2 {
        return 0;
    }
    let mut wn: i32 = 0;
    for i in 0..n {
        let (x0, y0) = polygon[i];
        let (x1, y1) = polygon[(i + 1) % n];
        if y0 <= py {
            if y1 > py && cross_sign((x0, y0), (x1, y1), (px, py)) > 0.0 {
                wn += 1;
            }
        } else if y1 <= py && cross_sign((x0, y0), (x1, y1), (px, py)) < 0.0 {
            wn -= 1;
        }
    }
    wn
}

/// Test whether two line segments properly intersect.
///
/// Uses the cross-product orientation test: segments (a1,a2) and (b1,b2)
/// intersect if and only if the endpoints of each segment lie on opposite
/// sides of the line containing the other segment.
pub(crate) fn segments_intersect(
    a1: (f64, f64),
    a2: (f64, f64),
    b1: (f64, f64),
    b2: (f64, f64),
) -> bool {
    let d1 = cross_sign(b1, b2, a1);
    let d2 = cross_sign(b1, b2, a2);
    let d3 = cross_sign(a1, a2, b1);
    let d4 = cross_sign(a1, a2, b2);
    d1 * d2 < 0.0 && d3 * d4 < 0.0
}

/// Compute the 2D cross product sign value `(b - a) × (c - a)`.
///
/// Returns a positive value if `a → b → c` is counterclockwise, negative
/// if clockwise, and zero if the points are collinear.
pub(crate) fn cross_sign(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}

#[cfg(test)]
#[cfg(test)]
#[path = "trim_tests.rs"]
mod tests;
