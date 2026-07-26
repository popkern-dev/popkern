//! Gates for the general n-sided fill primitive.
//!
//! The fixtures here are built from plain geometry rather than from a
//! fillet pipeline, which is the point: the primitive must work for a
//! caller that has only a ring of edges and the surfaces beside them.

use super::*;
use crate::geometry::surface::SurfaceEval;

/// A regular n-gon opening in the z = 0 plane, surrounded by a
/// notional flat neighbour whose normal is +z everywhere.
///
/// With a planar surround and a planar boundary, the transversal
/// directions all point inward and the fill should be a shallow,
/// well-behaved cap.
fn planar_opening(n: usize, radius: f64) -> Vec<FillBoundaryEdge> {
    let mut edges = Vec::with_capacity(n);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    for k in 0..n {
        let a0 = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        let a1 = 2.0 * std::f64::consts::PI * ((k as f64) + 0.5) / (n as f64);
        let a2 = 2.0 * std::f64::consts::PI * ((k + 1) as f64) / (n as f64);
        let p = |a: f64| Point3::new(radius * f64::cos(a), radius * f64::sin(a), 0.0);
        let (start, mid, end) = (p(a0), p(a1), p(a2));
        // Boundary tangent: along the loop, so the transversal
        // (normal x tangent) points into the opening.
        let tangent =
            |from: Point3, to: Point3| (to - from).normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
        let t = tangent(start, end);
        edges.push(FillBoundaryEdge {
            start,
            mid,
            end,
            normal_start: normal,
            normal_mid: normal,
            normal_end: normal,
            tangent_start: t,
            tangent_mid: t,
            tangent_end: t,
            extra_samples: Vec::new(),
        });
    }
    edges
}

#[test]
fn fills_every_supported_valence() {
    for n in 3..=12 {
        let boundary = planar_opening(n, 1.0);
        let fill = fill_n_sided_boundary(&boundary, 1e-9)
            .unwrap_or_else(|e| panic!("n={n} planar opening must fill, got {e}"));
        assert_eq!(fill.surfaces.len(), n, "one patch per side at n={n}");
        assert!(
            fill.worst_boundary_gap < 1e-12,
            "n={n} fill must interpolate its own boundary, worst gap {:.3e}",
            fill.worst_boundary_gap
        );
        for (k, surface) in fill.surfaces.iter().enumerate() {
            let p = surface.point_at(0.5, 0.5);
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "n={n} patch {k} evaluated to a non-finite point"
            );
        }
    }
}

#[test]
fn patches_meet_at_a_single_pole() {
    let fill = fill_n_sided_boundary(&planar_opening(5, 1.0), 1e-9).expect("pentagon fill");
    for (k, surface) in fill.surfaces.iter().enumerate() {
        let corner = surface.point_at(0.0, 0.0);
        assert!(
            corner.distance_to(&fill.pole) < 1e-9,
            "patch {k} corner {corner:?} does not sit at the pole {:?}",
            fill.pole
        );
    }
}

#[test]
fn too_few_and_too_many_sides_decline() {
    let two = planar_opening(3, 1.0)[..2].to_vec();
    assert!(matches!(
        fill_n_sided_boundary(&two, 1e-9),
        Err(NSidedFillError::UnsupportedSideCount { got: 2 })
    ));
    let thirteen = planar_opening(13, 1.0);
    assert!(matches!(
        fill_n_sided_boundary(&thirteen, 1e-9),
        Err(NSidedFillError::UnsupportedSideCount { got: 13 })
    ));
}

#[test]
fn open_boundary_declines_naming_the_gap() {
    let mut boundary = planar_opening(4, 1.0);
    // Pull one side's end away so the ring no longer closes.
    boundary[1].end = Point3::new(5.0, 5.0, 5.0);
    match fill_n_sided_boundary(&boundary, 1e-6) {
        Err(NSidedFillError::BoundaryNotClosed { side, gap, .. }) => {
            assert_eq!(side, 1, "the moved side must be the one reported");
            assert!(gap > 1.0, "the reported gap {gap:.3e} must be the real one");
        }
        other => panic!("an open boundary must decline, got {other:?}"),
    }
}

#[test]
fn degenerate_frame_declines_rather_than_folding() {
    let mut boundary = planar_opening(4, 1.0);
    // Make the normal parallel to the tangent: no transversal exists.
    boundary[2].normal_mid = boundary[2].tangent_mid;
    match fill_n_sided_boundary(&boundary, 1e-9) {
        Err(NSidedFillError::DegenerateFrame { side, sample, .. }) => {
            assert_eq!(side, 2);
            assert_eq!(sample, 1, "the mid sample is the one that was broken");
        }
        other => panic!("a degenerate frame must decline, got {other:?}"),
    }
}

#[test]
fn zero_length_side_declines() {
    // Collapse side 0 to a point while keeping the ring closed, so
    // the closure check (which runs first, by design) passes and the
    // degeneracy is what actually declines.
    let mut boundary = planar_opening(4, 1.0);
    let p = boundary[0].start;
    boundary[0].mid = p;
    boundary[0].end = p;
    boundary[1].start = p;
    assert!(matches!(
        fill_n_sided_boundary(&boundary, 1e-6),
        Err(NSidedFillError::DegenerateSide { side: 0, .. })
    ));
}

/// A non-planar opening — the case a real hole presents. The rim is
/// tilted so the surround's normals genuinely vary around the loop.
#[test]
fn non_planar_opening_fills_and_interpolates() {
    let n = 6;
    let radius = 1.0;
    let mut edges = Vec::with_capacity(n);
    let height = |a: f64| 0.3 * f64::sin(2.0 * a);
    for k in 0..n {
        let a0 = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        let a1 = 2.0 * std::f64::consts::PI * ((k as f64) + 0.5) / (n as f64);
        let a2 = 2.0 * std::f64::consts::PI * ((k + 1) as f64) / (n as f64);
        let p = |a: f64| Point3::new(radius * f64::cos(a), radius * f64::sin(a), height(a));
        let (start, mid, end) = (p(a0), p(a1), p(a2));
        let t = (end - start)
            .normalized()
            .expect("tilted rim side has length");
        // Surround normal: outward radial crossed with the tangent,
        // i.e. a plausible tangent plane for a wall around the hole.
        let radial = Vec3::new(f64::cos(a1), f64::sin(a1), 0.0);
        let normal = t
            .cross(&radial)
            .normalized()
            .expect("frame is non-degenerate");
        edges.push(FillBoundaryEdge {
            start,
            mid,
            end,
            normal_start: normal,
            normal_mid: normal,
            normal_end: normal,
            tangent_start: t,
            tangent_mid: t,
            tangent_end: t,
            extra_samples: Vec::new(),
        });
    }
    let fill = fill_n_sided_boundary(&edges, 1e-9).expect("tilted rim must fill");
    assert_eq!(fill.surfaces.len(), n);
    assert!(
        fill.worst_boundary_gap < 1e-12,
        "fill must interpolate a non-planar boundary too, worst gap {:.3e}",
        fill.worst_boundary_gap
    );
}
