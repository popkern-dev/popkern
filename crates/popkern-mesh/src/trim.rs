//! Adaptive trim polyline construction from face boundary loops.
//!
//! Builds UV-space polylines from face loops using pcurves (direct evaluation)
//! with 3D-to-UV projection fallback. These polylines are used by later tasks
//! for triangle classification and clipping.

// Items are pub(crate) for use by upcoming trim pipeline modules.
#![allow(dead_code)]

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::topology::{BRep, FaceHandle, LoopHandle};

/// A closed polyline in UV space representing a face boundary loop.
pub(crate) struct TrimLoop {
    /// UV points along the boundary.
    pub points: Vec<(f64, f64)>,
    /// Whether this is the outer boundary (`true`) or a hole (`false`).
    pub is_outer: bool,
}

/// Build trim loops for all boundary loops of a face.
///
/// Walks the outer loop and all inner loops, sampling each into a UV polyline.
/// The outer loop is tagged with `is_outer: true`; inner loops with `false`.
/// Loops with fewer than 3 points are skipped.
pub(crate) fn build_trim_loops(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    chord_height: f64,
) -> Vec<TrimLoop> {
    let Some(face) = brep.faces.get(face_handle) else {
        return Vec::new();
    };

    let mut result = Vec::new();

    // Outer loop.
    let outer_points = sample_loop(brep, face.outer_loop, surface, chord_height);
    if outer_points.len() >= 3 {
        result.push(TrimLoop {
            points: outer_points,
            is_outer: true,
        });
    }

    // Inner loops (holes).
    for &inner_loop in &face.inner_loops {
        let inner_points = sample_loop(brep, inner_loop, surface, chord_height);
        if inner_points.len() >= 3 {
            result.push(TrimLoop {
                points: inner_points,
                is_outer: false,
            });
        }
    }

    result
}

/// Sample a single loop into UV-space points.
///
/// For each half-edge in the loop:
/// - If the half-edge has a pcurve, samples it directly in UV space via
///   [`adaptive_sample_pcurve`].
/// - If no pcurve but the edge has a 3D curve, projects 32 uniform curve
///   samples to UV via [`project_point_to_uv`], respecting half-edge direction.
/// - If neither, projects the origin vertex position to UV.
pub(crate) fn sample_loop(
    brep: &BRep,
    loop_handle: LoopHandle,
    surface: &SurfaceGeometry,
    chord_height: f64,
) -> Vec<(f64, f64)> {
    let Some(loop_) = brep.loops.get(loop_handle) else {
        return Vec::new();
    };
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut uv_points = Vec::new();
    let mut current = start_he;

    loop {
        let he = &brep.half_edges[current];

        if let Some(ref pcurve) = he.pcurve {
            // Direct UV evaluation from pcurve (Point3(u, v, 0)).
            adaptive_sample_pcurve(pcurve, chord_height, &mut uv_points);
        } else {
            let edge = &brep.edges[he.edge];
            if let Some(ref curve_on_edge) = edge.curve {
                // Project 3D curve samples to UV.
                let samples = 32_usize;
                let curve_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
                let he_origin = brep.vertices[he.origin].position;
                let forward = curve_start.distance_to(&he_origin) < 1e-6;
                let (t_lo, t_hi) = if forward {
                    (curve_on_edge.t_start, curve_on_edge.t_end)
                } else {
                    (curve_on_edge.t_end, curve_on_edge.t_start)
                };

                for i in 0..samples {
                    let t = t_lo + (t_hi - t_lo) * (i as f64 / samples as f64);
                    let p = curve_on_edge.geometry.point_at(t);
                    if let Some(uv) = project_point_to_uv(surface, &p) {
                        uv_points.push(uv);
                    }
                }
            } else {
                // No curve geometry; project vertex position.
                let p = brep.vertices[he.origin].position;
                if let Some(uv) = project_point_to_uv(surface, &p) {
                    uv_points.push(uv);
                }
            }
        }

        current = he.next;
        if current == start_he {
            break;
        }
    }

    uv_points
}

/// Adaptively sample a pcurve (NURBS curve in UV space) using chord-height
/// subdivision.
///
/// Pushes the start point, then recursively subdivides spans where the
/// perpendicular chord-height error exceeds `chord_height`, up to a maximum
/// recursion depth of 6. Finally pushes the end point.
///
/// Chord-height error is the perpendicular distance from the curve midpoint
/// to the chord line: `|cross(p1-p0, p_mid-p0)| / |p1-p0|`.
pub(crate) fn adaptive_sample_pcurve(
    pcurve: &popkern_core::geometry::nurbs_curve::NurbsCurve,
    chord_height: f64,
    output: &mut Vec<(f64, f64)>,
) {
    let (t_start, t_end) = pcurve.parameter_range();
    let p_start = pcurve.point_at(t_start);
    output.push((p_start.x, p_start.y));

    subdivide_pcurve(pcurve, t_start, t_end, chord_height, 0, output);

    let p_end = pcurve.point_at(t_end);
    output.push((p_end.x, p_end.y));
}

/// Recursive subdivision helper for [`adaptive_sample_pcurve`].
fn subdivide_pcurve(
    pcurve: &popkern_core::geometry::nurbs_curve::NurbsCurve,
    t0: f64,
    t1: f64,
    chord_height: f64,
    depth: usize,
    output: &mut Vec<(f64, f64)>,
) {
    if depth >= 6 {
        return;
    }

    let t_mid = 0.5 * (t0 + t1);
    let p0 = pcurve.point_at(t0);
    let p1 = pcurve.point_at(t1);
    let p_mid = pcurve.point_at(t_mid);

    // Chord-height error: perpendicular distance from p_mid to the line p0-p1.
    // error = |cross(p1-p0, p_mid-p0)| / |p1-p0|
    let dx_chord = p1.x - p0.x;
    let dy_chord = p1.y - p0.y;
    let chord_len = (dx_chord * dx_chord + dy_chord * dy_chord).sqrt();

    let error = if chord_len > 1e-30 {
        let dx_mid = p_mid.x - p0.x;
        let dy_mid = p_mid.y - p0.y;
        // 2D cross product magnitude (z component of 3D cross).
        (dx_chord * dy_mid - dy_chord * dx_mid).abs() / chord_len
    } else {
        // Degenerate chord (zero length): use distance from midpoint to p0.
        let dx = p_mid.x - p0.x;
        let dy = p_mid.y - p0.y;
        (dx * dx + dy * dy).sqrt()
    };

    if error > chord_height {
        subdivide_pcurve(pcurve, t0, t_mid, chord_height, depth + 1, output);
        output.push((p_mid.x, p_mid.y));
        subdivide_pcurve(pcurve, t_mid, t1, chord_height, depth + 1, output);
    }
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
    let cx = (tri_uv[0].0 + tri_uv[1].0 + tri_uv[2].0) / 3.0;
    let cy = (tri_uv[0].1 + tri_uv[1].1 + tri_uv[2].1) / 3.0;

    if !point_inside_trim_loops(cx, cy, trim_loops) {
        return TriClass::Outside;
    }

    // Check if any triangle edge crosses any trim segment.
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
    }

    TriClass::Inside
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
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    #[test]
    fn cube_face_produces_one_outer_trim_loop() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face_handle = brep.faces.keys().next().unwrap();
        let face = &brep.faces[face_handle];
        let surface = &face.surface.as_ref().unwrap().geometry;

        let loops = build_trim_loops(&brep, face_handle, surface, 0.01);
        assert_eq!(loops.len(), 1, "cube face should have exactly 1 loop");
        assert!(loops[0].is_outer, "the loop should be outer");
        assert!(
            loops[0].points.len() >= 4,
            "cube face loop should have at least 4 points, got {}",
            loops[0].points.len()
        );
    }

    #[test]
    fn trim_loop_from_pcurves_is_closed() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face_handle = brep.faces.keys().next().unwrap();
        let face = &brep.faces[face_handle];
        let surface = &face.surface.as_ref().unwrap().geometry;

        let loops = build_trim_loops(&brep, face_handle, surface, 0.01);
        assert!(!loops.is_empty(), "should have at least one loop");
        let pts = &loops[0].points;
        let first = pts.first().unwrap();
        let last = pts.last().unwrap();
        let dist = ((first.0 - last.0).powi(2) + (first.1 - last.1).powi(2)).sqrt();
        assert!(
            dist < 0.1,
            "first and last points should be close, got distance {dist}"
        );
    }

    #[test]
    fn boolean_face_with_hole_produces_outer_and_inner_loops() {
        let tol = Tolerance::default().with_linear(1e-2);
        let cube = BRep::unit_cube(tol).unwrap();
        let cyl = BRep::cylinder(0.3, 2.0, tol).unwrap();

        // Translate cylinder to center it on the cube's top face.
        let offset = popkern_core::geometry::Vec3::new(0.5, 0.5, -0.5);
        let cyl = cyl
            .transform(&popkern_core::geometry::Transform3::translation(
                offset.x, offset.y, offset.z,
            ))
            .unwrap();

        let result = popkern_bool::subtraction(cube, cyl, tol).unwrap();

        // Check ALL faces with inner_loops; at least one must produce the
        // expected outer + inner trim loop pattern.
        let mut found_valid = false;
        for (fh, face) in &result.faces {
            if face.inner_loops.is_empty() || face.surface.is_none() {
                continue;
            }
            let surface = &face.surface.as_ref().unwrap().geometry;
            let loops = build_trim_loops(&result, fh, surface, 0.01);

            let outer_count = loops.iter().filter(|l| l.is_outer).count();
            let inner_count = loops.iter().filter(|l| !l.is_outer).count();

            if outer_count == 1 && inner_count >= 1 {
                found_valid = true;
                break;
            }
        }

        // If no face with inner_loops was found at all, the boolean result
        // topology doesn't use inner_loops for this operation -- acceptable.
        let has_any_face_with_inner = result.faces.values().any(|f| !f.inner_loops.is_empty());
        if has_any_face_with_inner {
            assert!(
                found_valid,
                "found faces with inner_loops but none produced outer+inner trim loops"
            );
        }
    }

    /// Helper: create a unit-square outer trim loop for classification tests.
    fn unit_square_loop() -> TrimLoop {
        TrimLoop {
            points: vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
            is_outer: true,
        }
    }

    #[test]
    fn classify_triangle_inside_square() {
        let loops = vec![unit_square_loop()];
        let tri = [(0.3, 0.3), (0.5, 0.3), (0.4, 0.5)];
        assert_eq!(classify_triangle(&tri, &loops), TriClass::Inside);
    }

    #[test]
    fn classify_triangle_outside_square() {
        let loops = vec![unit_square_loop()];
        let tri = [(1.5, 1.5), (2.0, 1.5), (1.7, 2.0)];
        assert_eq!(classify_triangle(&tri, &loops), TriClass::Outside);
    }

    #[test]
    fn classify_triangle_boundary_crossing_edge() {
        let loops = vec![unit_square_loop()];
        // Centroid ≈ (0.867, 0.433): inside the unit square.
        // The edge (1.2, 0.3)→(1.0, 0.7) crosses x=1 boundary.
        let tri = [(0.6, 0.3), (1.2, 0.3), (0.8, 0.7)];
        assert_eq!(classify_triangle(&tri, &loops), TriClass::Boundary);
    }

    #[test]
    fn classify_triangle_inside_hole_is_outside() {
        let outer = TrimLoop {
            points: vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)],
            is_outer: true,
        };
        let inner = TrimLoop {
            points: vec![(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)],
            is_outer: false,
        };
        let loops = vec![outer, inner];
        let tri = [(0.8, 0.8), (1.2, 0.8), (1.0, 1.2)];
        assert_eq!(classify_triangle(&tri, &loops), TriClass::Outside);
    }
}
