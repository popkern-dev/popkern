//! PK-380 — honest diagnostics for an oversized batch fillet.
//!
//! r = 0.9 on all 12 edges of a unit cube passes every per-edge radius
//! gate (a single 0.9 fillet fits a 1.0-wide face), the crossed strips
//! of each opposite pair are ridge-completed (PK-159), and the corner
//! construction then fails. Before this task the payload said
//! "24 of 12 edges failed" (both trim sides of each edge reported) and
//! classified every edge `InvalidTopology` → "heal the input BRep" on
//! a perfectly valid cube. The honest payload is one report per edge,
//! `RadiusTooLarge { requested: 0.9, max_radius: Some(0.5) }` — the
//! proportional-shrink bound r·gap/(r_a+r_b) — with `ReduceRadius`.
use popkern_blend::{
    BlendError, BlendFailure, FilletOptions, FilletRadius, FilletSpec, RemediationAction,
    fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn all_edge_specs(cube: &BRep, radius: f64) -> Vec<FilletSpec> {
    cube.edges
        .keys()
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect()
}

#[test]
fn oversized_cube_fillet_reports_radius_too_large_once_per_edge() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let specs = all_edge_specs(&cube, 0.9);
    let Err(err) = fillet_edges(&cube, &specs, &FilletOptions::default(), tol, None) else {
        panic!("r=0.9 on every edge of a unit cube cannot construct");
    };
    let BlendError::PartialFailure {
        failed,
        total,
        failures,
        ..
    } = &err
    else {
        panic!("expected PartialFailure, got {err}");
    };
    assert_eq!((*failed, *total), (12, 12), "{err}");
    assert_eq!(failures.len(), 12, "one report per fillet edge");
    let mut edges: Vec<_> = failures.iter().map(|(e, _)| *e).collect();
    edges.sort();
    edges.dedup();
    assert_eq!(edges.len(), 12, "no edge reported twice");
    for (edge, failure) in failures {
        let BlendFailure::RadiusTooLarge {
            requested,
            max_radius,
            location,
        } = failure
        else {
            panic!("{edge:?}: expected RadiusTooLarge, got {failure:?}");
        };
        assert_eq!(*requested, 0.9);
        let max = max_radius.expect("bound is known");
        assert!((max - 0.5).abs() < 1e-12, "{edge:?}: bound {max} != gap/2");
        let p = location.expect("located at the input edge midpoint");
        // Unit-cube edge midpoints have exactly two coordinates in
        // {0, 1} and one equal to 0.5.
        let coords = [p.x, p.y, p.z];
        assert_eq!(coords.iter().filter(|c| **c == 0.5).count(), 1, "{p:?}");
        assert_eq!(
            coords.iter().filter(|c| **c == 0.0 || **c == 1.0).count(),
            2,
            "{p:?}"
        );
        let hint = failure.remediation_hint();
        assert_eq!(hint.action, RemediationAction::ReduceRadius);
        let suggested = hint.value.expect("concrete radius");
        assert!(suggested < 0.5 && suggested > 0.45, "suggested {suggested}");
    }
    assert_eq!(err.remediation_hint().action, RemediationAction::Precheck);
}

#[test]
fn partial_failure_constructor_enforces_count_invariant() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let edge = cube.edges.keys().next().expect("edge");
    let dup = |reason: &str| {
        (
            edge,
            BlendFailure::InvalidTopology {
                reason: reason.to_string(),
                location: None,
            },
        )
    };
    let err = BlendError::partial_failure(3, cube.clone(), vec![dup("first"), dup("second")]);
    let BlendError::PartialFailure {
        failed,
        total,
        failures,
        ..
    } = err
    else {
        panic!("constructor builds PartialFailure");
    };
    assert_eq!((failed, total), (1, 3));
    assert!(
        matches!(&failures[0].1, BlendFailure::InvalidTopology { reason, .. } if reason == "first")
    );
}

/// The ridge-completion class itself still serves: two opposite
/// edges of a 1.0-wide slab at r = 0.6 cross mid-face and complete by
/// the parallel-cylinder seam (PK-159). Attribution must only touch
/// failures.
#[test]
fn ridge_completed_pair_still_serves() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    // Two parallel top edges: y = 0 and y = 1 at z = 1.
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .filter(|&e| {
            let Some(c) = cube.edges.get(e).and_then(|edge| edge.curve.as_ref()) else {
                return false;
            };
            use popkern_core::geometry::curve::CurveEval;
            let a = c.geometry.point_at(c.t_start);
            let b = c.geometry.point_at(c.t_end);
            a.z == 1.0 && b.z == 1.0 && a.y == b.y
        })
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.6),
        })
        .collect();
    assert_eq!(specs.len(), 2);
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol, None)
        .expect("PK-159 ridge completion serves 0.6 + 0.6 across a 1.0 face");
    assert!(result.brep.validate().is_ok());
}
