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
//! Every gate here is numeric: validation (structure, check 9, open
//! loops, self-intersection) plus the exact volume, from a closed form
//! for orthogonal fixtures and from an independent half-space clip for
//! convex ones (a planar chamfer of a convex body is its intersection
//! with the bevel half-spaces). Scale and rotation sweeps show the
//! constructions are not fixture-specific.

#[path = "support/mesh_gauges.rs"]
#[allow(dead_code)] // the welded-mesh gauge is not used here
mod mesh_gauges;

use popkern_blend::{
    BlendError, BlendFailure, ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges,
};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

const D: f64 = 0.1;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// A prism over the CCW outline `xy` from `z0` up `h`.
fn prism(xy: &[(f64, f64)], z0: f64, h: f64) -> BRep {
    let outline: Vec<Point3> = xy.iter().map(|&(x, y)| p(x, y, z0)).collect();
    BRep::prism(&outline, Vec3::new(0.0, 0.0, 1.0), h, Tolerance::default()).expect("prism")
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
}

fn subtract(body: BRep, tool: BRep) -> BRep {
    popkern_bool::subtraction(body, tool, Tolerance::default(), None).expect("fixture subtract")
}

/// The 2 × 2 × 1 block with a 1 × 1 × 0.5 pocket in its top.
fn pocket_block() -> BRep {
    subtract(
        prism(&rect(-1.0, -1.0, 1.0, 1.0), 0.0, 1.0),
        prism(&rect(-0.5, -0.5, 0.5, 0.5), 0.5, 1.0),
    )
}

/// The edge with endpoints `a` and `b` (either order).
fn edge(brep: &BRep, a: Point3, b: Point3) -> EdgeHandle {
    let gate = 1e-9 * (1.0 + a.distance_to(&b));
    brep.edges
        .iter()
        .find(|(_, e)| {
            let he = &brep.half_edges[e.half_edge];
            let (o, d) = (
                brep.vertices[he.origin].position,
                brep.vertices[brep.half_edges[he.twin].origin].position,
            );
            (o.distance_to(&a) < gate && d.distance_to(&b) < gate)
                || (o.distance_to(&b) < gate && d.distance_to(&a) < gate)
        })
        .map(|(h, _)| h)
        .unwrap_or_else(|| panic!("no edge {a:?} – {b:?}"))
}

fn chamfer(brep: &BRep, edges: &[(Point3, Point3)], d: f64) -> Result<BRep, BlendError> {
    let specs: Vec<ChamferSpec> = edges
        .iter()
        .map(|&(a, b)| ChamferSpec {
            edge: edge(brep, a, b),
            kind: ChamferKind::Equal(d),
        })
        .collect();
    chamfer_edges(
        brep,
        &specs,
        &ChamferOptions::default(),
        Tolerance::default(),
        None,
    )
    .map(|r| r.brep)
}

/// Valid by every validator and exactly `expected` in volume.
fn assert_exact(label: &str, brep: &BRep, expected: f64) {
    brep.validate()
        .unwrap_or_else(|e| panic!("{label}: invalid: {e}"));
    assert_eq!(
        mesh_gauges::brep_open_defects(brep),
        (0, 0, 0),
        "{label}: validate_all defects"
    );
    let volume = mesh_gauges::mesh_volume(brep);
    let rel = (volume - expected).abs() / expected;
    assert!(
        rel < 1e-9,
        "{label}: volume {volume:.12} vs exact {expected:.12} (rel {rel:.1e})"
    );
}

/// Removed volume of an equal-distance chamfer ring around a
/// rectangular opening `w × h` (four reflex corners): the frustum of
/// the inset rectangles minus the opening's own prism.
fn opening_ring(w: f64, h: f64, d: f64) -> f64 {
    (w + h) * d * d + 4.0 * d * d * d / 3.0
}

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

/// The reflex pocket-rim pair and the non-concurrent hexagon pair at
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
    }
}

/// A concave (fill) chamfer and a convex (cut) chamfer ending at the
/// same vertex have no construction: the decline names the corner and
/// points at chamfering them separately, and hands the input back.
#[test]
fn mixed_corner_declines_naming_the_corner() {
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
    let err = chamfer(
        &l,
        &[
            (p(1.0, 1.0, 0.0), p(1.0, 1.0, 0.5)),
            (p(2.0, 1.0, 0.5), p(1.0, 1.0, 0.5)),
        ],
        D,
    )
    .expect_err("mixed corner declines");
    let BlendError::PartialFailure {
        failures,
        partial_result,
        ..
    } = err
    else {
        panic!("expected a per-edge decline, got {err:?}");
    };
    assert_eq!(
        partial_result.faces.len(),
        l.faces.len(),
        "input handed back"
    );
    assert!(!failures.is_empty());
    for (_, failure) in &failures {
        let BlendFailure::VertexBlendFailed {
            reason, location, ..
        } = failure
        else {
            panic!("expected the corner to be named, got {failure:?}");
        };
        assert!(reason.starts_with("mixed chamfer corner"), "{reason}");
        assert!(location.is_some_and(|q| q.distance_to(&p(1.0, 1.0, 0.5)) < 1e-9));
        let hint = failure.remediation_hint();
        assert!(
            hint.rationale.contains("separate operations"),
            "{}",
            hint.rationale
        );
    }
}
