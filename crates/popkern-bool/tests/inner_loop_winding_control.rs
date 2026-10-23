//! PK-923 negative control: a drilled plate whose bore wall AND both
//! hole loops are reversed together. Every twin still pairs and every
//! face's outer loop still agrees with its flag, so the structural
//! checks and the winding-honesty audit pass, but the bore's material
//! sits on the wrong side: the class PK-916's half-inverted union
//! served. `validate_all` must name both hole loops.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::error::ValidationError;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

fn loop_half_edges(b: &BRep, l: LoopHandle) -> Vec<HalfEdgeHandle> {
    let start = b.loops[l].half_edge;
    let mut out = vec![start];
    let mut h = b.half_edges[start].next;
    while h != start {
        out.push(h);
        h = b.half_edges[h].next;
    }
    out
}

/// Reverse `loops` (whose edges both lie in `loops`) and flip the
/// flags of `faces`, as `BRep::reverse_body_orientation` does for a
/// whole body.
fn reverse_region(b: &mut BRep, loops: &[LoopHandle], faces: &[FaceHandle]) {
    let hes: Vec<HalfEdgeHandle> = loops.iter().flat_map(|&l| loop_half_edges(b, l)).collect();
    let origins: Vec<_> = hes
        .iter()
        .map(|&h| b.half_edges[b.half_edges[h].twin].origin)
        .collect();
    for (&h, o) in hes.iter().zip(origins) {
        let he = &mut b.half_edges[h];
        he.origin = o;
        std::mem::swap(&mut he.next, &mut he.prev);
        he.pcurve = he.pcurve.as_ref().map(|p| p.reversed());
    }
    for &h in &hes {
        let v = b.half_edges[h].origin;
        b.vertices[v].outgoing = h;
        let e = b.half_edges[h].edge;
        if b.edges[e].half_edge == b.half_edges[h].twin {
            continue;
        }
        if b.edges[e].half_edge == h && hes.contains(&b.half_edges[h].twin) {
            b.edges[e].half_edge = b.half_edges[h].twin;
        }
    }
    for &f in faces {
        if let Some(s) = b.faces[f].surface.as_mut() {
            s.orientation = !s.orientation;
        }
    }
}

#[test]
fn a_hole_wound_with_its_outer_loop_is_rejected() {
    let tol = Tolerance::default();
    let square = [
        Point3::new(-2.0, -2.0, 0.0),
        Point3::new(2.0, -2.0, 0.0),
        Point3::new(2.0, 2.0, 0.0),
        Point3::new(-2.0, 2.0, 0.0),
    ];
    let mut plate = BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 1.0, tol).expect("plate");
    plate.reverse_body_orientation();
    let drill = BRep::cylinder(0.5, 3.0, tol)
        .expect("drill")
        .transform(&popkern_core::geometry::Transform3::translation(
            0.0, 0.0, -1.0,
        ))
        .expect("drill placed");
    let body = boolean(BoolOperation::Subtraction, plate, drill, tol, None).expect("drilled");
    assert!(
        body.validate_all().is_empty(),
        "positive control: {:?}",
        body.validate_all()
    );

    // The bore wall: the faces on the cylinder; the hole loops: the
    // caps' inner loops.
    let wall: Vec<FaceHandle> = body
        .faces
        .iter()
        .filter(|(_, f)| {
            f.surface.as_ref().is_some_and(|s| {
                matches!(
                    s.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }
                )
            })
        })
        .map(|(h, _)| h)
        .collect();
    let holes: Vec<LoopHandle> = body
        .faces
        .values()
        .flat_map(|f| f.inner_loops.clone())
        .collect();
    assert!(
        !wall.is_empty() && holes.len() == 2,
        "wall {wall:?} holes {holes:?}"
    );
    let mut bad = body.clone();
    let mut loops = holes.clone();
    for &f in &wall {
        loops.push(bad.faces[f].outer_loop);
        loops.extend(bad.faces[f].inner_loops.clone());
    }
    reverse_region(&mut bad, &loops, &wall);

    assert!(
        bad.validate().is_ok(),
        "structural checks pass: {:?}",
        bad.validate().err()
    );
    let errors = bad.validate_all();
    let flagged: Vec<LoopHandle> = errors
        .iter()
        .filter_map(|e| match e {
            ValidationError::InnerLoopWinding { loop_, .. } => Some(*loop_),
            _ => None,
        })
        .collect();
    let mut want = holes.clone();
    want.sort();
    let mut got = flagged.clone();
    got.sort();
    assert_eq!(got, want, "all errors: {errors:?}");
    assert!(
        !errors
            .iter()
            .any(|e| matches!(e, ValidationError::WindingDishonest { .. })),
        "the outer loops stay honest: {errors:?}"
    );
    // The boolean's own result gate declines it too, so a serve of
    // this class falls to the rescue instead of shipping.
    let gate = popkern_bool::validate::post_check::post_check(&bad, &tol, false);
    assert!(
        gate.findings
            .iter()
            .any(|f| matches!(f.error, Some(ValidationError::InnerLoopWinding { .. }))),
        "post_check passes the inverted bore: {:?}",
        gate.findings
    );
    assert!(!gate.is_ok());
    assert!(popkern_bool::validate::post_check::post_check(&body, &tol, false).is_ok());
}

/// PK-924: a seamless cylinder band, its top rim an inner ring of the
/// lateral face (the seam killed by `kemr`, the form a STEP file may
/// carry). The ring wraps the period, so it has no UV winding and the
/// inner-loop check skips it, and the band validates clean. A side
/// witness by even-odd trim containment was tried and FALSIFIED on
/// exactly this body: ray parity across a seamless band reads its open
/// UV polylines and called the honest band's interior "outside". The
/// unmeasured class is a registry row (decline-class-registry.md,
/// "Inner-loop sense on loops without UV winding"); this pin keeps any
/// future witness honest on the band.
#[test]
fn a_seamless_band_validates_clean() {
    let tol = Tolerance::default();
    let mut band = BRep::cylinder(1.0, 1.0, tol).expect("cylinder");
    let seam = band
        .half_edges
        .iter()
        .find(|(_, h)| band.half_edges[h.twin].loop_ == h.loop_)
        .map(|(k, _)| k)
        .expect("seam half-edge");
    let ring = band.kemr(seam).expect("kemr").inner_loop;
    let lateral = band.loops[ring].face;
    let s = band.faces[lateral]
        .surface
        .as_ref()
        .expect("surface")
        .geometry
        .clone();
    assert!(
        band.loop_signed_uv_area(ring, &s).is_err(),
        "the ring has no UV winding"
    );
    assert!(band.validate_all().is_empty(), "{:?}", band.validate_all());
}
