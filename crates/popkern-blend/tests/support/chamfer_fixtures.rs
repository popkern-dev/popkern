//! Shared fixtures and exact gates for the chamfer corner tests
//! (PK-994, PK-996). Include beside `support/mesh_gauges.rs`.

use popkern_blend::{BlendError, ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

pub(crate) const D: f64 = 0.1;

pub(crate) fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// A prism over the CCW outline `xy` from `z0` up `h`.
pub(crate) fn prism(xy: &[(f64, f64)], z0: f64, h: f64) -> BRep {
    let outline: Vec<Point3> = xy.iter().map(|&(x, y)| p(x, y, z0)).collect();
    BRep::prism(&outline, Vec3::new(0.0, 0.0, 1.0), h, Tolerance::default()).expect("prism")
}

pub(crate) fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
}

pub(crate) fn subtract(body: BRep, tool: BRep) -> BRep {
    popkern_bool::subtraction(body, tool, Tolerance::default(), None).expect("fixture subtract")
}

/// The 2 × 2 × 1 block with a 1 × 1 × 0.5 pocket in its top.
pub(crate) fn pocket_block() -> BRep {
    subtract(
        prism(&rect(-1.0, -1.0, 1.0, 1.0), 0.0, 1.0),
        prism(&rect(-0.5, -0.5, 0.5, 0.5), 0.5, 1.0),
    )
}

/// The edge with endpoints `a` and `b` (either order).
pub(crate) fn edge(brep: &BRep, a: Point3, b: Point3) -> EdgeHandle {
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

pub(crate) fn chamfer(brep: &BRep, edges: &[(Point3, Point3)], d: f64) -> Result<BRep, BlendError> {
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
pub(crate) fn assert_exact(label: &str, brep: &BRep, expected: f64) {
    brep.validate()
        .unwrap_or_else(|e| panic!("{label}: invalid: {e}"));
    assert_eq!(
        super::mesh_gauges::brep_open_defects(brep),
        (0, 0, 0),
        "{label}: validate_all defects"
    );
    let volume = super::mesh_gauges::mesh_volume(brep);
    let rel = (volume - expected).abs() / expected;
    assert!(
        rel < 1e-9,
        "{label}: volume {volume:.12} vs exact {expected:.12} (rel {rel:.1e})"
    );
}

/// Removed volume of an equal-distance chamfer ring around a
/// rectangular opening `w × h` (four reflex corners): the frustum of
/// the inset rectangles minus the opening's own prism.
pub(crate) fn opening_ring(w: f64, h: f64, d: f64) -> f64 {
    (w + h) * d * d + 4.0 * d * d * d / 3.0
}

pub(crate) fn l_prism() -> BRep {
    prism(
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
    )
}

/// Fills first (PK-996): a concave chamfer of height `h` adds `h·d²/2`;
/// a convex chamfer of length `len` ending against the fill removes its
/// wedge minus the part the fill's plane keeps, `len·d²/2 − 2d³/3`.
pub(crate) fn fill(h: f64) -> f64 {
    h * D * D / 2.0
}
pub(crate) fn cut_against_fill(len: f64) -> f64 {
    len * D * D / 2.0 - 2.0 * D * D * D / 3.0
}
