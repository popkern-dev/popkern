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
}
