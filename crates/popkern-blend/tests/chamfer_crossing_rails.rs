//! PK-1006 / PK-1007: bevels whose rails cross past a face's side.
//!
//! A face whose rails leave it no region of the pipeline's shape (a
//! cone's slant triangle chamfered on all three sides past its
//! inradius) is consumed: the request must be served by bevel growth,
//! exactly, never by cutting the crossing rails pairwise (that served
//! a validator-clean body with material outside a bevel plane). Wide
//! bevels on prisms that reach the multi-chord trim's spur re-home
//! must serve exactly or decline with the input back, never panic.
//! Every serve gates on validation and the half-space clip volume.

#[path = "support/apex_cones.rs"]
#[allow(dead_code)] // not every cone is used here
mod apex_cones;
#[path = "support/chamfer_fixtures.rs"]
#[allow(dead_code)] // not every fixture is used here
mod chamfer_fixtures;
#[path = "support/halfspace_oracle.rs"]
mod halfspace_oracle;
#[path = "support/mesh_gauges.rs"]
#[allow(dead_code)] // the welded-mesh gauge is not used here
mod mesh_gauges;

use apex_cones::*;
use chamfer_fixtures::*;
use popkern_blend::{BlendError, ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Chamfer `edges` by `d`: `Some(result)` when served, `None` when the
/// request declined with the input back.
fn served_or_declined(
    body: &BRep,
    edges: &[(Point3, Point3)],
    d: f64,
    label: &str,
) -> Option<BRep> {
    let specs: Vec<ChamferSpec> = edges
        .iter()
        .map(|&(a, b)| ChamferSpec {
            edge: edge(body, a, b),
            kind: ChamferKind::Equal(d),
        })
        .collect();
    match chamfer_edges(
        body,
        &specs,
        &ChamferOptions::default(),
        Tolerance::default(),
        None,
    ) {
        Ok(result) => Some(result.brep),
        Err(BlendError::PartialFailure { partial_result, .. }) => {
            assert_eq!(
                partial_result.faces.len(),
                body.faces.len(),
                "{label}: the input comes back"
            );
            None
        }
        Err(e) => panic!("{label}: unexpected error {e}"),
    }
}

/// Every slant triangle of the 5-sided cone, chamfered on its three
/// sides 15% past its inradius, is consumed and serves exactly.
#[test]
fn slant_triangle_beveled_past_its_inradius_is_consumed() {
    let body = cone5();
    let (apex, ends) = apex_edges(&body);
    for i in 0..ends.len() {
        let (b, c) = (ends[i], ends[(i + 1) % ends.len()]);
        let area = 0.5 * (b - apex).cross(&(c - apex)).length();
        let half_perimeter =
            0.5 * (apex.distance_to(&b) + b.distance_to(&c) + c.distance_to(&apex));
        let d = 1.15 * area / half_perimeter;
        let edges = [(apex, b), (b, c), (c, apex)];
        let label = format!("5-cone face {i}, d = {d:.4}");
        let served = served_or_declined(&body, &edges, d, &label)
            .unwrap_or_else(|| panic!("{label}: declined"));
        assert_exact(
            &label,
            &served,
            halfspace_oracle::chamfered_volume(&body, &edges, d),
        );
    }
}

/// Wide bevels on a quadrilateral and a triangular prism, found by a
/// randomized audit: their rails cross past a side, and the multi-chord
/// trim indexed an anchor an earlier trim's spur re-home had removed
/// (a panic). The consumption gate now declines them before any trim
/// and growth serves the quadrilateral: exact or a clean decline.
#[test]
fn wide_prism_bevels_serve_exactly_or_decline() {
    // Outlines in the audit's generation order (from the corner
    // nearest azimuth 0), which fixes the topology's handle order.
    let quad = [
        (0.8339831305420121, -0.1513835059876938),
        (0.01108700165153909, 0.8475387343998855),
        (-0.8475976167483774, -0.004807077724734521),
        (-0.21746518142731328, -0.8192397224137552),
    ];
    let tri = [
        (1.4214380136237683, 0.3368727704002522),
        (-0.6980401805476629, 1.28324167495948),
        (-1.0825173896201936, -0.9808799066051392),
    ];
    let mut served = 0;
    for (name, outline, h, d, picks) in [
        (
            "quad",
            &quad[..],
            0.9213001226518827,
            0.47172784305286775,
            // (from corner, to corner, top?) per edge; `None` = upright.
            vec![
                (3, 0, Some(false)),
                (1, 2, Some(true)),
                (1, 2, Some(false)),
                (1, 1, None),
                (2, 3, Some(true)),
            ],
        ),
        (
            "triangle",
            &tri[..],
            0.7924981214153894,
            0.43719125768871386,
            vec![
                (1, 2, Some(false)),
                (0, 1, Some(true)),
                (2, 0, Some(true)),
                (0, 0, None),
                (2, 0, Some(false)),
                (1, 2, Some(true)),
            ],
        ),
    ] {
        let body = prism(outline, 0.0, h);
        let at = |k: usize, z: f64| p(outline[k].0, outline[k].1, z);
        let edges: Vec<(Point3, Point3)> = picks
            .iter()
            .map(|&(i, j, top)| match top {
                Some(top) => {
                    let z = if top { h } else { 0.0 };
                    (at(i, z), at(j, z))
                }
                None => (at(i, 0.0), at(i, h)),
            })
            .collect();
        if let Some(result) = served_or_declined(&body, &edges, d, name) {
            assert_exact(
                name,
                &result,
                halfspace_oracle::chamfered_volume(&body, &edges, d),
            );
            served += 1;
        }
    }
    assert!(served >= 1, "the quadrilateral serves through growth");
}
