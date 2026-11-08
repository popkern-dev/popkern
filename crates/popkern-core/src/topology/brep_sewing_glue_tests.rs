use super::super::brep_sewing::{BrepSewingProvenance, sew_clean_boundary_edge_pairs};
use super::*;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BodyKind;
use popkern_math::PortableMath;

/// One square sheet on the hinge x in [0, 1], leaving it toward
/// (0, cos θ, sin θ): a lamina on a wire square, its loop running the
/// hinge along +x when `forward`. Returns the sheet and its hinge edge.
fn page(theta: f64, forward: bool) -> (BRep, EdgeHandle) {
    let mut b = BRep::new(Tolerance::default());
    let out = b.outside;
    let w = Vec3::new(0.0, theta.pcos(), theta.psin());
    let x = Vec3::new(1.0, 0.0, 0.0);
    let o = Point3::new(0.0, 0.0, 0.0);
    let corners = if forward {
        [o, o + x, o + x + w, o + w]
    } else {
        [o, o + w, o + x + w, o + x]
    };
    let (u_axis, v_axis) = if forward { (x, w) } else { (w, x) };
    let v: Vec<_> = corners
        .iter()
        .map(|&p| b.make_acorn(out, p).expect("acorn").vertex)
        .collect();
    let wires: Vec<_> = (0..4)
        .map(|i| b.make_wire_edge(out, v[i], v[(i + 1) % 4]).expect("wire"))
        .collect();
    let cycle: Vec<_> = wires.iter().map(|r| r.fin).collect();
    let face = b
        .make_face_from_wire_cycle(
            &cycle,
            Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Plane {
                    origin: o,
                    normal: u_axis.cross(&v_axis),
                    u_axis,
                    v_axis,
                },
                orientation: true,
            }),
        )
        .expect("face")
        .face;
    b.recompute_face_pcurves(face).expect("pcurves");
    (b, wires[if forward { 0 } else { 3 }].edge)
}

fn options() -> BrepSewingOptions {
    BrepSewingOptions {
        tolerance: 1e-6,
        samples: 3,
        require_surface_agreement: false,
        validate_output: true,
        allow_nonmanifold_edges: true,
    }
}

fn provenance() -> BrepSewingProvenance {
    BrepSewingProvenance {
        first: "page".to_string(),
        second: "page".to_string(),
    }
}

/// Three sheets on one hinge (an odd fan: two run it one way, one the
/// other): sewn pairwise into one edge, the third onto the edge the
/// first pair already joined, a ring of three in radial order; one
/// general body that validates.
#[test]
fn three_sheets_sew_onto_one_hinge_edge() {
    let mut body = BRep::new(Tolerance::default());
    let mut hinges = Vec::new();
    for (theta, forward) in [
        (0.0, true),
        (120.0_f64.to_radians(), false),
        (240.0_f64.to_radians(), true),
    ] {
        let (p, hinge) = page(theta, forward);
        let map = body.append_body(&p);
        hinges.push(map.edges[&hinge]);
    }
    let pairs = [
        BrepSewingEdgePair {
            first: hinges[0],
            second: hinges[1],
            provenance: provenance(),
        },
        BrepSewingEdgePair {
            first: hinges[0],
            second: hinges[2],
            provenance: provenance(),
        },
    ];
    let sewn = sew_clean_boundary_edge_pairs(&body, &pairs, options()).expect("sewn");
    assert_eq!(
        sewn.report.accepted_pairs, 2,
        "{:?}",
        sewn.report.declined_pairs
    );
    let b = sewn.brep;
    assert_eq!(b.radial_fins(hinges[0]).len(), 3);
    assert_eq!(b.body_kind(), BodyKind::General);
    assert_eq!(
        b.radial_order(hinges[0], &b.tolerance).expect("order"),
        b.radial_fins(hinges[0])
    );
    let errors = b.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
}

/// Without the non-manifold option the third sheet declines as before.
#[test]
fn manifold_sewing_still_declines_the_third_sheet() {
    let mut body = BRep::new(Tolerance::default());
    let mut hinges = Vec::new();
    for (theta, forward) in [
        (0.0, true),
        (120.0_f64.to_radians(), false),
        (240.0_f64.to_radians(), true),
    ] {
        let (p, hinge) = page(theta, forward);
        let map = body.append_body(&p);
        hinges.push(map.edges[&hinge]);
    }
    let pairs = [
        BrepSewingEdgePair {
            first: hinges[0],
            second: hinges[1],
            provenance: provenance(),
        },
        BrepSewingEdgePair {
            first: hinges[0],
            second: hinges[2],
            provenance: provenance(),
        },
    ];
    let manifold = BrepSewingOptions {
        allow_nonmanifold_edges: false,
        validate_output: false,
        ..options()
    };
    let sewn = sew_clean_boundary_edge_pairs(&body, &pairs, manifold).expect("sewn");
    assert_eq!(sewn.report.accepted_pairs, 1);
    assert_eq!(sewn.report.declined_pairs.len(), 1);
}
