//! Torus trims stay on both surfaces end to end (PK-941).
//!
//! The junction refine pinned each clipped record end to the nearest
//! face-cell iso-boundary and moved only its end control point. A face
//! spanning its whole period has no cell, so a cylinder-seam end was
//! dragged 1.1e-2 along the trim onto a torus patch boundary, and the
//! curve's end span bowed 1.1e-4 off both surfaces (3.6e-4 for a
//! sphere bite) while its end point and the validator's samples sat on
//! the trim. Every periodic chart's seam is now a pin candidate.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::geometry::analytic_surfaces::analytic_signed_distance;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn rotations() -> Vec<(&'static str, Transform3)> {
    vec![
        ("unrotated", Transform3::identity()),
        ("z 0.7", Transform3::rotation_z(0.7)),
        ("x 1.1", Transform3::rotation_x(1.1)),
        ("y -0.4", Transform3::rotation_y(-0.4)),
        (
            "mixed",
            Transform3::rotation_x(0.3)
                .compose(&Transform3::rotation_y(1.9))
                .compose(&Transform3::rotation_z(-2.2)),
        ),
    ]
}

/// The farthest any edge between two different surfaces strays from
/// either, over 4000 samples of its curve.
fn worst_trim_deviation(body: &BRep) -> f64 {
    let mut worst: f64 = 0.0;
    for edge in body.edges.values() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let h = &body.half_edges[edge.half_edge];
        let surfaces: Vec<_> = [h.loop_, body.half_edges[h.twin].loop_]
            .iter()
            .filter_map(|&l| body.faces[body.loops[l].face].surface.as_ref())
            .map(|s| &s.geometry)
            .collect();
        if surfaces.len() != 2 || surfaces[0].kind_name() == surfaces[1].kind_name() {
            continue;
        }
        for k in 0..4000 {
            let f = (k as f64 + 0.5) / 4000.0;
            let p = curve
                .geometry
                .point_at(curve.t_start + (curve.t_end - curve.t_start) * f);
            for s in &surfaces {
                let d = analytic_signed_distance(s, &p).expect("analytic").abs();
                worst = worst.max(d);
            }
        }
    }
    worst
}

#[test]
fn drilled_and_bitten_tori_keep_their_trims_on_both_surfaces() {
    let tol = Tolerance::default();
    let tools = [
        (
            "drilled through the tube",
            BRep::cylinder(0.15, 2.0, tol)
                .and_then(|c| c.transform(&Transform3::translation(1.0, 0.0, -1.0))),
        ),
        (
            "bitten at the outer equator",
            BRep::sphere(0.2, tol)
                .and_then(|s| s.transform(&Transform3::translation(1.3, 0.0, 0.0))),
        ),
    ];
    for (case, tool) in tools {
        let tool = tool.expect("tool");
        for (name, r) in rotations() {
            let torus = BRep::torus(1.0, 0.3, tol)
                .and_then(|t| t.transform(&r))
                .expect("torus");
            let cut = boolean(
                BoolOperation::Subtraction,
                torus,
                tool.transform(&r).expect("rotate"),
                tol,
                None,
            )
            .expect("boolean");
            let worst = worst_trim_deviation(&cut);
            assert!(
                worst < 10.0 * tol.linear(),
                "{case}, {name}: a trim strays {worst:e} off its surfaces"
            );
        }
    }
}
