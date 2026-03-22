//! General surface-curve intersection fallback.
//!
//! Orchestrates: overlap detection -> curve Bezier decomposition -> surface
//! bounding-box grid (8x8) -> AABB culling -> Newton 3x3 polish -> dedup -> sort.

use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::intersection::surface_curve_newton::newton_surface_curve;
use crate::geometry::intersection::surface_curve_overlap::detect_surface_curve_overlaps;
use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::Point3;
use crate::tolerance::Tolerance;

/// Number of cells per axis in the surface bounding-box grid.
const GRID_N: usize = 8;

// ---------------------------------------------------------------------------
// Surface AABB grid
// ---------------------------------------------------------------------------

struct SurfaceCell {
    u_range: (f64, f64),
    v_range: (f64, f64),
    min: Point3,
    max: Point3,
}

/// Build an `n x n` grid of surface cells, each with an AABB computed from
/// the four corner evaluations.
fn build_surface_grid(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    n: usize,
) -> Vec<SurfaceCell> {
    let (u0, u1) = u_range;
    let (v0, v1) = v_range;
    let du = (u1 - u0) / n as f64;
    let dv = (v1 - v0) / n as f64;

    // Pre-evaluate all grid corner points: (n+1) x (n+1).
    let mut pts = vec![Point3::new(0.0, 0.0, 0.0); (n + 1) * (n + 1)];
    for iu in 0..=n {
        for iv in 0..=n {
            let u = u0 + iu as f64 * du;
            let v = v0 + iv as f64 * dv;
            pts[iu * (n + 1) + iv] = surface.point_at(u, v);
        }
    }

    let mut cells = Vec::with_capacity(n * n);
    for iu in 0..n {
        for iv in 0..n {
            let corners = [
                pts[iu * (n + 1) + iv],
                pts[iu * (n + 1) + iv + 1],
                pts[(iu + 1) * (n + 1) + iv],
                pts[(iu + 1) * (n + 1) + iv + 1],
            ];
            let mut bmin = corners[0];
            let mut bmax = corners[0];
            for c in &corners[1..] {
                if c.x < bmin.x { bmin.x = c.x; }
                if c.y < bmin.y { bmin.y = c.y; }
                if c.z < bmin.z { bmin.z = c.z; }
                if c.x > bmax.x { bmax.x = c.x; }
                if c.y > bmax.y { bmax.y = c.y; }
                if c.z > bmax.z { bmax.z = c.z; }
            }
            cells.push(SurfaceCell {
                u_range: (u0 + iu as f64 * du, u0 + (iu + 1) as f64 * du),
                v_range: (v0 + iv as f64 * dv, v0 + (iv + 1) as f64 * dv),
                min: bmin,
                max: bmax,
            });
        }
    }
    cells
}

// ---------------------------------------------------------------------------
// Curve AABB
// ---------------------------------------------------------------------------

struct CurveSegment {
    t_range: (f64, f64),
    min: Point3,
    max: Point3,
}

/// Compute AABB for a NURBS curve by sampling a set of points.
fn curve_aabb(nurbs: &NurbsCurve, t_range: (f64, f64)) -> (Point3, Point3) {
    let n_samples = 16;
    let (t0, t1) = t_range;
    let mut p = nurbs.point_at(t0);
    let mut bmin = p;
    let mut bmax = p;
    for i in 1..=n_samples {
        let frac = i as f64 / n_samples as f64;
        p = nurbs.point_at(t0 + frac * (t1 - t0));
        if p.x < bmin.x { bmin.x = p.x; }
        if p.y < bmin.y { bmin.y = p.y; }
        if p.z < bmin.z { bmin.z = p.z; }
        if p.x > bmax.x { bmax.x = p.x; }
        if p.y > bmax.y { bmax.y = p.y; }
        if p.z > bmax.z { bmax.z = p.z; }
    }
    (bmin, bmax)
}

/// Check AABB overlap with a tolerance margin.
fn aabb_overlap(
    a_min: &Point3, a_max: &Point3,
    b_min: &Point3, b_max: &Point3,
    margin: f64,
) -> bool {
    a_min.x - margin <= b_max.x && a_max.x + margin >= b_min.x
        && a_min.y - margin <= b_max.y && a_max.y + margin >= b_min.y
        && a_min.z - margin <= b_max.z && a_max.z + margin >= b_min.z
}

// ---------------------------------------------------------------------------
// Bezier decomposition
// ---------------------------------------------------------------------------

/// Decompose a NURBS curve into Bezier segments by splitting at all interior
/// knots (knots with multiplicity < degree are raised to full multiplicity,
/// then the result is split).
fn decompose_to_bezier_segments(nurbs: &NurbsCurve) -> Vec<(NurbsCurve, (f64, f64))> {
    let (t_min, t_max) = nurbs.parameter_range();
    let degree = nurbs.degree;

    // Collect unique interior knot values.
    let mut interior_knots: Vec<f64> = Vec::new();
    for &k in &nurbs.knots {
        if k > t_min && k < t_max
            && interior_knots.last().is_none_or(|&last| (k - last).abs() > 1e-14) {
            interior_knots.push(k);
        }
    }

    if interior_knots.is_empty() {
        // Already a single Bezier segment.
        return vec![(nurbs.clone(), (t_min, t_max))];
    }

    // For each interior knot, count its current multiplicity and insert enough
    // times to raise it to `degree`.  We collect all insertions.
    let mut knots_to_insert: Vec<f64> = Vec::new();
    for &k in &interior_knots {
        let mult = nurbs.knots.iter().filter(|&&v| (v - k).abs() < 1e-14).count();
        let needed = degree.saturating_sub(mult);
        for _ in 0..needed {
            knots_to_insert.push(k);
        }
    }

    let refined = if knots_to_insert.is_empty() {
        nurbs.clone()
    } else {
        match nurbs.insert_knots(&knots_to_insert) {
            Ok(c) => c,
            Err(_) => nurbs.clone(),
        }
    };

    // Now split at each interior knot.
    let mut segments: Vec<(NurbsCurve, (f64, f64))> = Vec::new();
    let mut remaining = refined;

    for &k in &interior_knots {
        let (r_min, _r_max) = remaining.parameter_range();
        match remaining.split_at(k) {
            Ok((left, right)) => {
                let left_range = left.parameter_range();
                segments.push((left, (r_min.max(left_range.0), k.min(left_range.1))));
                remaining = right;
            }
            Err(_) => {
                // If split fails (e.g. knot at boundary), skip.
            }
        }
    }

    // Add the last segment.
    let rem_range = remaining.parameter_range();
    segments.push((remaining, rem_range));

    segments
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// General surface-curve intersection solver.
///
/// Orchestrates: overlap detection, curve Bezier decomposition, surface
/// bounding-box grid, AABB culling, Newton 3x3 polish, dedup, sort.
pub(crate) fn intersect_surface_curve_general(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    curve: &CurveGeometry,
    t_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceCurveIntersection> {
    // --- Step 1: overlap detection ---
    let mut results = detect_surface_curve_overlaps(
        surface, u_range, v_range, curve, t_range, tolerance,
    );

    // --- Step 2: promote curve to NURBS and decompose into Bezier segments ---
    let nurbs = curve.to_nurbs_in_range(t_range);
    let bezier_segments = decompose_to_bezier_segments(&nurbs);

    // Compute AABB for each Bezier segment.
    let curve_segments: Vec<CurveSegment> = bezier_segments
        .into_iter()
        .map(|(seg, tr)| {
            let (bmin, bmax) = curve_aabb(&seg, tr);
            CurveSegment { t_range: tr, min: bmin, max: bmax }
        })
        .collect();

    // --- Step 3: build surface bbox grid ---
    let surface_cells = build_surface_grid(surface, u_range, v_range, GRID_N);

    // --- Step 4: AABB culling + Newton polish ---
    let margin = tolerance.linear * 10.0; // small margin for AABB overlap test
    let mut point_results: Vec<(f64, f64, f64, Point3)> = Vec::new();

    for cseg in &curve_segments {
        for cell in &surface_cells {
            if !aabb_overlap(&cseg.min, &cseg.max, &cell.min, &cell.max, margin) {
                continue;
            }

            // Seed Newton from the midpoints of the curve segment and surface cell.
            let t_mid = 0.5 * (cseg.t_range.0 + cseg.t_range.1);
            let u_mid = 0.5 * (cell.u_range.0 + cell.u_range.1);
            let v_mid = 0.5 * (cell.v_range.0 + cell.v_range.1);

            if let Some((t, u, v, pt)) =
                newton_surface_curve(curve, t_mid, surface, u_mid, v_mid, tolerance)
            {
                // Check that the result is within the requested ranges.
                if t >= t_range.0 - tolerance.parametric
                    && t <= t_range.1 + tolerance.parametric
                    && u >= u_range.0 - tolerance.parametric
                    && u <= u_range.1 + tolerance.parametric
                    && v >= v_range.0 - tolerance.parametric
                    && v <= v_range.1 + tolerance.parametric
                {
                    point_results.push((t, u, v, pt));
                }
            }
        }
    }

    // --- Step 5: deduplicate ---
    let dedup_dist = tolerance.linear * 10.0;
    let mut unique: Vec<(f64, f64, f64, Point3)> = Vec::new();
    'outer: for (t, u, v, pt) in point_results {
        for &(et, _eu, _ev, ref ep) in &unique {
            if (t - et).abs() < dedup_dist || pt.distance_to(ep) < dedup_dist {
                continue 'outer;
            }
        }
        unique.push((t, u, v, pt));
    }

    // Also skip points that fall inside any overlap range.
    for (t, u, v, pt) in unique {
        let in_overlap = results.iter().any(|r| {
            if let SurfaceCurveIntersection::Overlap { t_range: (ta, tb) } = r {
                t >= *ta && t <= *tb
            } else {
                false
            }
        });
        if !in_overlap {
            results.push(SurfaceCurveIntersection::Point { point: pt, t, u, v });
        }
    }

    // --- Step 6: sort by t ---
    results.sort_by(|a, b| {
        let ta = match a {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range } => t_range.0,
        };
        let tb = match b {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            SurfaceCurveIntersection::Overlap { t_range } => t_range.0,
        };
        ta.partial_cmp(&tb).unwrap_or(std::cmp::Ordering::Equal)
    });

    results
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::CurveGeometry;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;

    // -----------------------------------------------------------------------
    // Test 1: NURBS line piercing a plane (verify against analytic line-plane)
    //
    // Plane z=0, line along z-axis from (0,0,-1) to (0,0,1).
    // Should find one intersection at t=0.5 (midpoint), point (0,0,0).
    // -----------------------------------------------------------------------
    #[test]
    fn nurbs_line_piercing_plane() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };

        // Build a NURBS curve that represents a line from (0,0,-1) to (0,0,1).
        let nurbs_line = CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, -1.0),
                Point3::new(0.0, 0.0, 1.0),
            ],
            weights: vec![1.0, 1.0],
        });

        let tol = Tolerance::default();
        let results = intersect_surface_curve_general(
            &plane,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &nurbs_line,
            (0.0, 1.0),
            &tol,
        );

        let points: Vec<_> = results
            .iter()
            .filter(|r| matches!(r, SurfaceCurveIntersection::Point { .. }))
            .collect();

        assert_eq!(points.len(), 1, "Expected exactly 1 point intersection, got {}", points.len());
        if let SurfaceCurveIntersection::Point { point, t, .. } = &points[0] {
            assert!(
                (point.z).abs() < 1e-5,
                "Point should be on z=0 plane, got z={}",
                point.z
            );
            assert!(
                (*t - 0.5).abs() < 1e-5,
                "t should be ~0.5, got {}",
                t
            );
        }
    }

    // -----------------------------------------------------------------------
    // Test 2: Circle (as NURBS) piercing a plane — expect 2 intersections.
    //
    // Plane z=0.  Circle in xz-plane centered at origin, radius 1.
    // The circle crosses z=0 at (1,0,0) and (-1,0,0).
    // -----------------------------------------------------------------------
    #[test]
    fn circle_piercing_plane_two_points() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };

        // Circle in the xz-plane: center=(0,0,0), normal=(0,1,0), radius=1.
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            radius: 1.0,
        };

        let tol = Tolerance::default();
        let results = intersect_surface_curve_general(
            &plane,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &circle,
            (0.0, std::f64::consts::TAU),
            &tol,
        );

        let points: Vec<_> = results
            .iter()
            .filter(|r| matches!(r, SurfaceCurveIntersection::Point { .. }))
            .collect();

        assert!(
            points.len() >= 2,
            "Expected at least 2 point intersections for circle piercing plane, got {}",
            points.len()
        );

        // All intersection points should be on the z=0 plane.
        for p in &points {
            if let SurfaceCurveIntersection::Point { point, .. } = p {
                assert!(
                    point.z.abs() < 1e-4,
                    "Intersection point should be on z=0 plane, got z={}",
                    point.z
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // Test 3: Non-intersecting curve and surface (empty).
    //
    // Plane z=0.  Line from (0,0,1) to (0,0,2) — entirely above the plane.
    // -----------------------------------------------------------------------
    #[test]
    fn non_intersecting_empty() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };

        let nurbs_line = CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 1.0),
                Point3::new(0.0, 0.0, 2.0),
            ],
            weights: vec![1.0, 1.0],
        });

        let tol = Tolerance::default();
        let results = intersect_surface_curve_general(
            &plane,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &nurbs_line,
            (0.0, 1.0),
            &tol,
        );

        let points: Vec<_> = results
            .iter()
            .filter(|r| matches!(r, SurfaceCurveIntersection::Point { .. }))
            .collect();

        assert!(
            points.is_empty(),
            "Expected no intersections, got {}",
            points.len()
        );
    }

    // -----------------------------------------------------------------------
    // Test 4: Invariant — curve.point_at(t) ~ surface.point_at(u, v)
    //
    // For every Point intersection returned, verify the invariant.
    // -----------------------------------------------------------------------
    #[test]
    fn invariant_curve_surface_match() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };

        let nurbs_line = CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, -1.0),
                Point3::new(0.0, 0.0, 1.0),
            ],
            weights: vec![1.0, 1.0],
        });

        let tol = Tolerance::default();
        let results = intersect_surface_curve_general(
            &plane,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &nurbs_line,
            (0.0, 1.0),
            &tol,
        );

        for r in &results {
            if let SurfaceCurveIntersection::Point { point, t, u, v } = r {
                let curve_pt = nurbs_line.point_at(*t);
                let surface_pt = plane.point_at(*u, *v);
                assert!(
                    point.distance_to(&curve_pt) < 1e-4,
                    "curve.point_at(t) should match intersection point"
                );
                assert!(
                    point.distance_to(&surface_pt) < 1e-4,
                    "surface.point_at(u,v) should match intersection point"
                );
            }
        }
    }
}
