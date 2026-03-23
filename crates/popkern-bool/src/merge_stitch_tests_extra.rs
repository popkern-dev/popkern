//! Extra stitch.rs coverage tests split out to keep file sizes under 500 lines.
#[cfg(test)]
use std::collections::HashSet;

use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

use crate::merge::merge_faces;
use crate::select::SelectionResult;
use crate::split::{FragmentMap, SplitResult};

fn identity_split(a: BRep, b: BRep) -> SplitResult {
    let fragments_a = FragmentMap {
        origin: a.faces.keys().map(|fh| (fh, fh)).collect(),
        intersection_edges: HashSet::new(),
    };
    let fragments_b = FragmentMap {
        origin: b.faces.keys().map(|fh| (fh, fh)).collect(),
        intersection_edges: HashSet::new(),
    };
    SplitResult {
        a,
        b,
        fragments_a,
        fragments_b,
    }
}

fn select_a_only(split: &SplitResult) -> SelectionResult {
    let keep_a: HashSet<_> = split.a.faces.keys().collect();
    SelectionResult {
        keep_a,
        keep_b: HashSet::new(),
        flip_a: HashSet::new(),
        flip_b: HashSet::new(),
    }
}

// ── stitch.rs coverage: refit_curve_to_vertices degenerate / non-line ─────────

/// Assign a zero-direction Line curve (degenerate) to a cube edge with
/// mismatched parameters. refit_curve_to_vertices returns None at line 334-336.
#[test]
fn stitch_refit_curve_degenerate_line_does_not_panic() {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};

    let tol = Tolerance::default();
    let mut cube = BRep::unit_cube(tol).unwrap();

    let eh = cube.edges.keys().next().unwrap();
    cube.edges[eh].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(0.0, 0.0, 0.0), // zero-length — degenerate
        },
        t_start: 100.0,
        t_end: 200.0,
    });

    let empty = BRep::new(tol);
    let split = identity_split(cube, empty);
    let selection = select_a_only(&split);
    let result = merge_faces(&split, &selection, &tol).unwrap();
    assert_eq!(result.faces.len(), 6);
}

/// Assign a Circle curve (non-Line) to a cube edge with mismatched parameters.
/// refit_curve_to_vertices hits the `_ => None` branch at line 355.
#[test]
fn stitch_refit_curve_non_line_does_not_panic() {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};

    let tol = Tolerance::default();
    let mut cube = BRep::unit_cube(tol).unwrap();

    let eh = cube.edges.keys().next().unwrap();
    cube.edges[eh].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        },
        t_start: 100.0,
        t_end: 200.0,
    });

    let empty = BRep::new(tol);
    let split = identity_split(cube, empty);
    let selection = select_a_only(&split);
    let result = merge_faces(&split, &selection, &tol).unwrap();
    assert_eq!(result.faces.len(), 6);
}
