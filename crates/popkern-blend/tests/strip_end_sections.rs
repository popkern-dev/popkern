//! Strip ends on oblique faces and at pyramid apexes (PK-530).
//!
//! A single fillet on a pyramid edge runs out at its two end vertices.
//! A base edge ends on the neighbouring slant faces, which are oblique
//! to its axis, so the end section is an ellipse. A slant edge ends at
//! the apex, into the fan of the other apex faces: one arc per face.
//! Before PK-530 the chain-end settlement assumed a single
//! perpendicular cap, declined everything else, and `fillet_edges`
//! returned `Ok` with the stale L of straight stubs through the
//! original corner: up to 8.4e-2 off the strip, caught only by check
//! 9 (hexagonal pyramid base edge). Now:
//! - cylinder strips settle exactly on any planar end (ellipse
//!   section), at any fan valence (split corner);
//! - fitted G2 strips settle on marched sections;
//! - anything still inconsistent declines typed
//!   (`InconsistentGeometry`), with the input handed back.

use popkern_blend::{
    BlendError, BlendFailure, ContinuityLevel, FilletOptions, FilletRadius, FilletSpec,
    fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

fn pyramid(n: usize) -> BRep {
    BRep::n_gonal_pyramid(n, 1.0, 1.2, Tolerance::default()).expect("pyramid")
}

/// A base edge (both ends at z = 0) or a slant edge (to the apex).
fn pick(b: &BRep, slant: bool) -> EdgeHandle {
    b.edges
        .iter()
        .find(|(_, e)| {
            let he = &b.half_edges[e.half_edge];
            let z0 = b.vertices[he.origin].position.z;
            let z1 = b.vertices[b.half_edges[he.twin].origin].position.z;
            if slant {
                (z0 - z1).abs() > 0.5
            } else {
                z0.abs() < 1e-9 && z1.abs() < 1e-9
            }
        })
        .map(|(h, _)| h)
        .expect("edge")
}

fn fillet(
    b: &BRep,
    edge: EdgeHandle,
    r: f64,
    continuity: ContinuityLevel,
) -> Result<popkern_blend::BlendResult, BlendError> {
    let options = FilletOptions {
        continuity,
        ..FilletOptions::default()
    };
    fillet_edges(
        b,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        }],
        &options,
        Tolerance::default(),
        None,
    )
}

fn assert_serves_clean(label: &str, result: Result<popkern_blend::BlendResult, BlendError>) {
    let result = result.unwrap_or_else(|e| panic!("{label}: must serve: {e}"));
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "{label}: validate_all: {errors:?}");
    assert!(
        result.brep.edges.values().all(|e| e.curve.is_some()),
        "{label}: every edge carries a curve"
    );
}

/// Cylinder strips: every base end (oblique cap) and every slant end
/// (the base cap and the apex fan, valence 3 to 6) settles exactly.
#[test]
fn single_edge_cylinder_strips_end_exactly_on_oblique_caps_and_apex_fans() {
    for n in 3..=6 {
        let b = pyramid(n);
        for slant in [false, true] {
            for r in [0.05, 0.12] {
                let label = format!("n={n} slant={slant} r={r}");
                assert_serves_clean(&label, fillet(&b, pick(&b, slant), r, ContinuityLevel::G1));
            }
        }
    }
}

/// Fitted G2 strips settle their oblique and fan ends on marched
/// sections: slant edges at every valence, and base edges where the
/// strip meets its rails (dihedral ≥ 59.5°).
#[test]
fn single_edge_g2_strips_end_on_marched_sections() {
    for n in 3..=6 {
        let b = pyramid(n);
        assert_serves_clean(
            &format!("slant n={n}"),
            fillet(&b, pick(&b, true), 0.07, ContinuityLevel::G2),
        );
    }
    for n in [3, 4] {
        let b = pyramid(n);
        assert_serves_clean(
            &format!("base n={n}"),
            fillet(&b, pick(&b, false), 0.07, ContinuityLevel::G2),
        );
    }
}

/// PK-531 pin: a G2 strip on a pentagonal or hexagonal pyramid's base
/// edge (dihedral ≤ 56°) disagrees with its own base rail (2.8e-2 to
/// 7.8e-2 along its length). It once served that body `Ok`, and the
/// territory gate now declines it typed, with the input handed back.
#[test]
fn g2_base_strips_at_acute_dihedrals_decline_typed() {
    for n in [5, 6] {
        let b = pyramid(n);
        let Err(BlendError::PartialFailure {
            partial_result,
            failures,
            ..
        }) = fillet(&b, pick(&b, false), 0.07, ContinuityLevel::G2)
        else {
            panic!("n={n}: PK-531 serves this now — flip the pin to a serve");
        };
        assert_eq!(
            partial_result.faces.len(),
            b.faces.len(),
            "input handed back"
        );
        assert!(
            failures
                .iter()
                .any(|(_, f)| matches!(f, BlendFailure::InconsistentGeometry { .. })),
            "n={n}: typed territory decline: {failures:?}"
        );
    }
}

/// Blend on blend (PK-534): after a cube's top edges are rounded, a
/// fillet on each vertical edge ends against the two top strips, which
/// meet along their PK-152 seam. That is a fan end with curved faces
/// and a curved fan edge. It once served the stale L, (√2 − 1)·r =
/// 4.1e-2 off. It now settles on marched strip ∩ strip arcs, with the
/// seam trimmed at its crossing, for G1 verticals on a G1 or G2 top.
#[test]
fn blend_on_blend_vertical_ends_settle_against_the_top_strips() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let top: Vec<FilletSpec> = cube
        .edges
        .iter()
        .filter(|(_, e)| {
            let he = &cube.half_edges[e.half_edge];
            cube.vertices[he.origin].position.z > 0.99
                && cube.vertices[cube.half_edges[he.twin].origin].position.z > 0.99
        })
        .map(|(edge, _)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.15),
        })
        .collect();
    for top_continuity in [ContinuityLevel::G1, ContinuityLevel::G2] {
        let rounded = fillet_edges(
            &cube,
            &top,
            &FilletOptions {
                continuity: top_continuity,
                ..FilletOptions::default()
            },
            tol,
            None,
        )
        .expect("top edges round")
        .brep;
        let verticals: Vec<EdgeHandle> = rounded
            .edges
            .iter()
            .filter(|(_, e)| {
                e.curve.as_ref().is_some_and(|c| {
                    use popkern_core::geometry::CurveEval;
                    let (a, b) = (c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end));
                    (a.x - b.x).abs() < 1e-9 && (a.y - b.y).abs() < 1e-9 && (a.z - b.z).abs() > 0.5
                })
            })
            .map(|(h, _)| h)
            .collect();
        assert_eq!(verticals.len(), 4);
        let specs: Vec<FilletSpec> = verticals
            .iter()
            .map(|&edge| FilletSpec {
                edge,
                radius: FilletRadius::Constant(0.1),
            })
            .collect();
        let result = fillet_edges(&rounded, &specs, &FilletOptions::default(), tol, None);
        let result =
            result.unwrap_or_else(|e| panic!("{top_continuity:?} top: verticals serve: {e}"));
        assert_eq!(
            result.repair_report.terminal_arc_trim.repaired_corners, 4,
            "{top_continuity:?} top: every vertical's top end settles"
        );
        let errors = result.brep.validate_all();
        assert!(
            errors.is_empty(),
            "{top_continuity:?} top: validate_all: {errors:?}"
        );
    }
}
