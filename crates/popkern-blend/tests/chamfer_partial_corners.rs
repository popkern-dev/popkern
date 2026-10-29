//! PK-994: chamfers that end at partially chamfered corners.
//!
//! Two corner classes declined on every input before PK-994:
//!
//! - **Reflex ends.** A chamfered edge whose end vertex also carries a
//!   concave edge (a pocket or hole rim, a side tunnel's mouth, a
//!   re-entrant L corner). The support's rail end lies on the
//!   neighbour edge's EXTENSION, so the end face grows to the bevel.
//! - **Bevel pairs.** Two chamfered edges meeting where a third edge
//!   stays sharp (an L of box top edges, a top edge with its vertical).
//!   The bevels mitre; when the four planes do not concur (hexagon,
//!   slanted ramp) the corner needs an extra connector edge.
//!
//! Mixed corners (PK-996) are `chamfer_mixed_corners.rs`.
//!
//! Every gate here is numeric: validation (structure, check 9, open
//! loops, self-intersection) plus the exact volume, from a closed form
//! for orthogonal fixtures and from an independent half-space clip for
//! convex ones (a planar chamfer of a convex body is its intersection
//! with the bevel half-spaces). Scale and rotation sweeps show the
//! constructions are not fixture-specific.

#[path = "support/chamfer_fixtures.rs"]
#[allow(dead_code)] // not every fixture is used here
mod chamfer_fixtures;
#[path = "support/mesh_gauges.rs"]
#[allow(dead_code)] // the welded-mesh gauge is not used here
mod mesh_gauges;

use chamfer_fixtures::*;
use popkern_blend::{ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn pocket_rim_edge_ends_flush_with_the_pocket_walls() {
    let body = pocket_block();
    let faces_before = body.faces.len();
    let out = chamfer(&body, &[(p(-0.5, 0.5, 1.0), p(0.5, 0.5, 1.0))], D).expect("serves");
    assert_exact("pocket rim edge", &out, 3.5 - D * D / 2.0);
    // One bevel; each end's triangle joins its pocket wall.
    assert_eq!(
        out.faces.len(),
        faces_before + 1,
        "terminal triangles merged"
    );
}

#[test]
fn pocket_rim_edge_two_distance_reflex_ends() {
    let body = pocket_block();
    let top = body
        .faces
        .keys()
        .find(|&f| {
            body.face_outer_loop_positions(f)
                .iter()
                .all(|q| (q.z - 1.0).abs() < 1e-12)
        })
        .expect("top face");
    let spec = ChamferSpec {
        edge: edge(&body, p(-0.5, 0.5, 1.0), p(0.5, 0.5, 1.0)),
        kind: ChamferKind::TwoDistance {
            face_1: top,
            dist_1: 0.1,
            dist_2: 0.2,
        },
    };
    let out = chamfer_edges(
        &body,
        &[spec],
        &ChamferOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("serves")
    .brep;
    assert_exact("pocket rim two-distance", &out, 3.5 - 0.5 * 0.1 * 0.2);
}

#[test]
fn pocket_rim_adjacent_edges_mitre_outside_the_corner() {
    let out = chamfer(
        &pocket_block(),
        &[
            (p(-0.5, 0.5, 1.0), p(0.5, 0.5, 1.0)),
            (p(0.5, 0.5, 1.0), p(0.5, -0.5, 1.0)),
        ],
        D,
    )
    .expect("serves");
    // Two wedges plus the mitred piece outside the reflex corner (d³/3).
    assert_exact(
        "pocket rim pair",
        &out,
        3.5 - 2.0 * D * D / 2.0 - D * D * D / 3.0,
    );
}

#[test]
fn pocket_rim_closed_loop() {
    let rim = [
        (p(-0.5, 0.5, 1.0), p(0.5, 0.5, 1.0)),
        (p(0.5, 0.5, 1.0), p(0.5, -0.5, 1.0)),
        (p(0.5, -0.5, 1.0), p(-0.5, -0.5, 1.0)),
        (p(-0.5, -0.5, 1.0), p(-0.5, 0.5, 1.0)),
    ];
    let out = chamfer(&pocket_block(), &rim, D).expect("serves");
    assert_exact("pocket rim loop", &out, 3.5 - opening_ring(1.0, 1.0, D));
}

#[test]
fn side_tunnel_mouth_ring() {
    let body = subtract(
        prism(&rect(-1.0, -1.0, 1.0, 1.0), 0.0, 1.0),
        prism(&rect(-1.5, -0.4, 0.5, 0.4), 0.3, 0.4),
    );
    let (a, b, c, e) = (
        p(-1.0, -0.4, 0.3),
        p(-1.0, 0.4, 0.3),
        p(-1.0, 0.4, 0.7),
        p(-1.0, -0.4, 0.7),
    );
    let one = chamfer(&body, &[(a, b)], D).expect("one mouth edge serves");
    assert_exact("tunnel mouth edge", &one, 3.52 - 0.8 * D * D / 2.0);
    let ring = chamfer(&body, &[(a, b), (b, c), (c, e), (e, a)], D).expect("mouth ring serves");
    assert_exact("tunnel mouth ring", &ring, 3.52 - opening_ring(0.8, 0.4, D));
}

#[test]
fn reentrant_corner_top_edge() {
    let l = prism(
        &[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 2.0),
            (0.0, 2.0),
        ],
        0.0,
        0.5,
    );
    let out = chamfer(&l, &[(p(2.0, 1.0, 0.5), p(1.0, 1.0, 0.5))], D).expect("serves");
    assert_exact("re-entrant top edge", &out, 1.5 - D * D / 2.0);
}

#[test]
fn box_bevel_pairs_concurrent() {
    let block = prism(&rect(-1.0, -1.0, 1.0, 1.0), 0.0, 1.0);
    let top_y = (p(-1.0, 1.0, 1.0), p(1.0, 1.0, 1.0));
    let top_x = (p(1.0, 1.0, 1.0), p(1.0, -1.0, 1.0));
    let vertical = (p(1.0, 1.0, 0.0), p(1.0, 1.0, 1.0));
    // Pair overlap at an orthogonal corner: d³/3.
    let pair = chamfer(&block, &[top_y, top_x], D).expect("top L serves");
    assert_exact("box top L", &pair, 4.0 - (2.0 * D * D - D * D * D / 3.0));
    let with_vertical = chamfer(&block, &[top_y, vertical], D).expect("top + vertical serves");
    assert_exact(
        "box top + vertical",
        &with_vertical,
        4.0 - (D * D + D * D / 2.0 - D * D * D / 3.0),
    );
}

#[path = "support/halfspace_oracle.rs"]
mod halfspace_oracle;

#[test]
fn hexagon_top_edge_with_vertical_non_concurrent() {
    let s = 3f64.sqrt() / 2.0;
    let hex = [
        (1.0, 0.0),
        (0.5, s),
        (-0.5, s),
        (-1.0, 0.0),
        (-0.5, -s),
        (0.5, -s),
    ];
    let body = prism(&hex, 0.0, 1.0);
    let edges = [
        (p(1.0, 0.0, 1.0), p(0.5, s, 1.0)),
        (p(1.0, 0.0, 0.0), p(1.0, 0.0, 1.0)),
    ];
    let out = chamfer(&body, &edges, D).expect("serves");
    let exact = halfspace_oracle::chamfered_volume(&body, &edges, D);
    assert_exact("hexagon top + vertical", &out, exact);
}

#[test]
fn slanted_ramp_top_edges_non_concurrent() {
    // Trapezoid section in the XZ plane, swept 1 along −y: slanted sides.
    let section = [
        p(0.0, 0.0, 0.0),
        p(2.0, 0.0, 0.0),
        p(1.6, 0.0, 1.0),
        p(0.4, 0.0, 1.0),
    ];
    let ramp = BRep::prism(
        &section,
        Vec3::new(0.0, -1.0, 0.0),
        1.0,
        Tolerance::default(),
    )
    .expect("ramp");
    for edges in [
        [
            (p(0.4, 0.0, 1.0), p(1.6, 0.0, 1.0)),
            (p(1.6, 0.0, 1.0), p(1.6, -1.0, 1.0)),
        ],
        [
            (p(1.6, 0.0, 1.0), p(1.6, -1.0, 1.0)),
            (p(1.6, 0.0, 1.0), p(2.0, 0.0, 0.0)),
        ],
    ] {
        let out = chamfer(&ramp, &edges, D).expect("serves");
        let exact = halfspace_oracle::chamfered_volume(&ramp, &edges, D);
        assert_exact("ramp pair", &out, exact);
    }
}

/// The reflex pocket-rim pair, the non-concurrent hexagon pair and the
/// mixed L corner at
/// other scales and rotated off every axis: the same corners, so the
/// volume scales as s³ (and matches the clip of the placed body).
#[test]
fn reflex_pair_scale_and_rotation_sweep() {
    let unit = 3.5 - 2.0 * D * D / 2.0 - D * D * D / 3.0;
    let rotation =
        Transform3::rotation_around_axis(p(0.1, -0.2, 0.3), Vec3::new(1.0, 2.0, 3.0), 0.7)
            .expect("rotation");
    for s in [1e-2, 1.0, 1e2] {
        let place = rotation.compose(&Transform3::scale_uniform(s));
        let body = pocket_block().transform(&place).expect("placed fixture");
        let at = |q: Point3| place.apply_point(&q);
        let edges = [
            (at(p(-0.5, 0.5, 1.0)), at(p(0.5, 0.5, 1.0))),
            (at(p(0.5, 0.5, 1.0)), at(p(0.5, -0.5, 1.0))),
        ];
        let out = chamfer(&body, &edges, D * s).unwrap_or_else(|e| panic!("scale {s}: {e}"));
        assert_exact(&format!("pocket pair, scale {s}"), &out, unit * s * s * s);
        // The non-concurrent hexagon corner, placed the same way.
        let r = 3f64.sqrt() / 2.0;
        let hex = [
            (1.0, 0.0),
            (0.5, r),
            (-0.5, r),
            (-1.0, 0.0),
            (-0.5, -r),
            (0.5, -r),
        ];
        let body = prism(&hex, 0.0, 1.0)
            .transform(&place)
            .expect("placed hexagon");
        let edges = [
            (at(p(1.0, 0.0, 1.0)), at(p(0.5, r, 1.0))),
            (at(p(1.0, 0.0, 0.0)), at(p(1.0, 0.0, 1.0))),
        ];
        let out = chamfer(&body, &edges, D * s).unwrap_or_else(|e| panic!("hex scale {s}: {e}"));
        let exact = halfspace_oracle::chamfered_volume(&body, &edges, D * s);
        assert_exact(&format!("hexagon pair, scale {s}"), &out, exact);
        // The mixed (fills-first) corner, placed the same way.
        let body = l_prism().transform(&place).expect("placed L");
        let edges = [
            (at(p(1.0, 1.0, 0.0)), at(p(1.0, 1.0, 0.5))),
            (at(p(2.0, 1.0, 0.5)), at(p(1.0, 1.0, 0.5))),
        ];
        let out = chamfer(&body, &edges, D * s).unwrap_or_else(|e| panic!("L scale {s}: {e}"));
        let unit_l = 1.5 + fill(0.5) - cut_against_fill(1.0);
        assert_exact(&format!("L mixed, scale {s}"), &out, unit_l * s * s * s);
    }
}

/// A quadrilateral prism whose top corners are 93° and 96°: a bevel
/// pair's two rails cross the sharp top edge 0.14% and 0.55% of `d`
/// apart, so the corner is non-concurrent by less than 1e-3. Every
/// subset of its four edges serves exactly at three sizes (before
/// PK-1008 the second rail end snapped onto the first through the
/// loose anchor reuse, building the corner as concurrent with a mitre
/// 1.9e-4 off its plane: a decline at d = 0.1).
#[test]
fn nearly_orthogonal_bevel_pairs_stay_non_concurrent() {
    let q = [
        (-0.042953069041007465, -1.4393661841422616),
        (1.4140660281527422, -0.27209786147014325),
        (-0.18693179576301835, 1.427822286534676),
        (-1.4295605077637037, 0.17313790120790618),
    ];
    let h = 1.1882629623909862;
    let body = prism(&q, 0.0, h);
    let at = |k: usize, z: f64| p(q[k].0, q[k].1, z);
    let all = [
        (at(3, 0.0), at(3, h)),
        (at(1, h), at(2, h)),
        (at(0, 0.0), at(0, h)),
        (at(3, h), at(0, h)),
    ];
    for mask in 1..16usize {
        let edges: Vec<(Point3, Point3)> = (0..4)
            .filter(|i| mask & (1 << i) != 0)
            .map(|i| all[i])
            .collect();
        for d in [0.03, 0.1, 0.475] {
            let label = format!("quad edges {mask:04b}, d = {d}");
            let out = chamfer(&body, &edges, d).unwrap_or_else(|e| panic!("{label}: {e}"));
            assert_exact(
                &label,
                &out,
                halfspace_oracle::chamfered_volume(&body, &edges, d),
            );
        }
    }
}
