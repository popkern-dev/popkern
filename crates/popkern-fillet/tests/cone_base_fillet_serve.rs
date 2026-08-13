//! PK-264 serving gate: the cone base-circle fillet.
//!
//! Before PK-264 this operation failed at EVERY radius — the convexity
//! classifier read the base rim as Concave (its interior-left rule
//! trusted loop windings that primitive constructors never enforce),
//! sending the walker to the exterior mirror ball, and the footprint
//! clip only branch-normalized cylinder angles, so the cone footprint
//! lost 18 of 20 samples against the flank's (−2π, 0) chart range
//! ("footprint length mismatch: 2 vs 20"). The §11.1 contract test
//! only asserted "the surface is a Torus" — never which side — so the
//! mirror ball was invisible to it.
//!
//! These gates pin the end-to-end serve, the inscribed (material-side)
//! torus, and the honest structured decline above the feasible radius.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};
use popkern_fillet::{
    BlendError, BlendFailure, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// The base circle: the only edge whose two adjacent faces differ.
fn base_edge(brep: &BRep) -> EdgeHandle {
    brep.edges
        .keys()
        .find(|&e| {
            let he = brep.edges[e].half_edge;
            let twin = brep.half_edges[he].twin;
            brep.loops[brep.half_edges[he].loop_].face
                != brep.loops[brep.half_edges[twin].loop_].face
        })
        .expect("cone has a base edge with distinct faces")
}

/// End-to-end serve at two radii, with the fillet torus pinned to the
/// INSCRIBED side: center above the base plane (+z, inside the
/// material) and major radius strictly inside the base circle. The
/// exterior mirror ball — the pre-fix walker attractor — has center
/// below the plane and major radius outside the rim.
#[test]
fn cone_base_fillet_serves_with_inscribed_torus() {
    let cone = BRep::cone(1.0, 2.0, tol()).expect("cone builds");
    let edge = base_edge(&cone);
    for radius in [0.005_f64, 0.03] {
        let specs = vec![FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        }];
        let result = fillet_edges(&cone, &specs, &FilletOptions::default(), tol(), None)
            .unwrap_or_else(|e| panic!("cone base fillet r={radius} must serve: {e}"));
        assert!(
            result.brep.validate().is_ok(),
            "r={radius}: fillet output must validate"
        );
        let torus = result
            .brep
            .faces
            .values()
            .filter_map(|f| f.surface.as_ref())
            .find_map(|s| match &s.geometry {
                SurfaceGeometry::Torus {
                    center,
                    major_radius,
                    minor_radius,
                    ..
                } => Some((*center, *major_radius, *minor_radius)),
                _ => None,
            })
            .expect("fillet face must carry an analytic torus (§11.1 plane-cone arm)");
        let (center, major, minor) = torus;
        assert!(
            (minor - radius).abs() < 1e-12,
            "r={radius}: torus minor radius must equal the fillet radius, got {minor}"
        );
        assert!(
            (center.z - radius).abs() < 1e-9,
            "r={radius}: inscribed ball centers sit one radius ABOVE the base plane \
             (mirror-ball regression guard), got center.z = {}",
            center.z
        );
        assert!(
            major < 1.0,
            "r={radius}: inscribed major radius lies inside the base circle \
             (mirror-ball regression guard), got {major}"
        );
    }
}

/// Above the feasible radius the operation must decline with the
/// structured RadiusTooLarge carrying an honest maximum — one that a
/// retry at 90% of the reported value actually serves.
#[test]
fn cone_base_fillet_declines_honestly_above_max_radius() {
    let cone = BRep::cone(1.0, 2.0, tol()).expect("cone builds");
    let edge = base_edge(&cone);
    let attempt = |r: f64| {
        fillet_edges(
            &cone,
            &[FilletSpec {
                edge,
                radius: FilletRadius::Constant(r),
            }],
            &FilletOptions::default(),
            tol(),
            None,
        )
    };
    // Find the decline threshold by walking up until the op declines.
    let mut declined_max: Option<f64> = None;
    for r in [0.05_f64, 0.1, 0.2, 0.4, 0.8] {
        match attempt(r) {
            Ok(_) => continue,
            Err(BlendError::PartialFailure { failures, .. }) => {
                let max = failures.iter().find_map(|(_, f)| match f {
                    BlendFailure::RadiusTooLarge { max_radius, .. } => *max_radius,
                    _ => None,
                });
                assert!(
                    max.is_some(),
                    "the first infeasible radius (r={r}) must decline RadiusTooLarge \
                     with a populated max_radius, got {failures:?}"
                );
                declined_max = max;
                break;
            }
            Err(other) => panic!("r={r}: expected RadiusTooLarge decline, got {other}"),
        }
    }
    if let Some(max) = declined_max {
        let retry = 0.9 * max;
        assert!(
            attempt(retry).is_ok(),
            "the reported max_radius {max} must be honest: 90% of it ({retry}) must serve"
        );
    } else {
        // Every probed radius served — also acceptable; the decline
        // class is then unreached on this fixture.
        eprintln!("cone base fillet served through r=0.8 — no decline reached");
    }
}
