// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//! PK-376 pins: lamina orientation flags are derived from the measured
//! loop winding, both faces share one chart, and the region's front
//! face is the one facing the sketch normal — for a loop authored in
//! either direction, with or without a hole.

use crate::{SketchOnPlane, SketchPlane};
use popkern_constraint::Sketch;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::history::TrackingContext;
use popkern_core::naming::{OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn rectangle(sk: &mut Sketch, x0: f64, y0: f64, x1: f64, y1: f64, clockwise: bool) {
    let mut corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)];
    if clockwise {
        corners.reverse();
    }
    let pts: Vec<_> = corners.iter().map(|(x, y)| sk.add_point(*x, *y)).collect();
    for i in 0..4 {
        sk.add_line(pts[i], pts[(i + 1) % 4]);
    }
}

fn build(sk: Sketch) -> (crate::SketchBrepResult, TrackingContext) {
    let mut sop = SketchOnPlane {
        sketch: sk,
        plane: SketchPlane::xy(),
        element_names: Default::default(),
    };
    let mut tracking = TrackingContext::new(OperationId(1));
    sop.assign_element_names(tracking.operation);
    let result = sop
        .to_brep(Tolerance::default(), Some(&mut tracking))
        .expect("sketch builds");
    (result, tracking)
}

/// Newell normal of a face's outer loop.
fn loop_newell(brep: &BRep, face: FaceHandle) -> Vec3 {
    let pts = brep.face_outer_loop_positions(face);
    let mut n = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..pts.len() {
        let a = pts[i];
        let b = pts[(i + 1) % pts.len()];
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    n
}

/// Flag-folded outward normal of a planar face.
fn outward_normal(brep: &BRep, face: FaceHandle) -> Vec3 {
    let sof = brep.faces[face].surface.as_ref().expect("surface");
    match &sof.geometry {
        SurfaceGeometry::Plane { normal, .. } => {
            if sof.orientation {
                *normal
            } else {
                -*normal
            }
        }
        other => panic!("lamina face must be planar, got {other:?}"),
    }
}

fn assert_lamina_honest(brep: &BRep, front: FaceHandle) {
    // Winding audit (pcurve-based) and 3-D Newell agree with the flag
    // on every face, and the strict validator is clean.
    let errs = brep.validate_all();
    assert!(errs.is_empty(), "{errs:?}");
    for v in brep.audit_winding_honesty() {
        assert_eq!(v.honest, Some(true), "{v:?}");
    }
    for (fh, _) in &brep.faces {
        let newell = loop_newell(brep, fh);
        let outward = outward_normal(brep, fh);
        assert!(
            newell.dot(&outward) > 0.0,
            "{fh:?}: loop winds against its flag-folded normal"
        );
    }
    // The front face faces +Z; its twin faces −Z; they share one chart.
    assert!(outward_normal(brep, front).z > 0.0);
    let back = brep.faces.keys().find(|f| *f != front).expect("back face");
    assert!(outward_normal(brep, back).z < 0.0);
    let geom = |f: FaceHandle| format!("{:?}", brep.faces[f].surface.as_ref().map(|s| &s.geometry));
    assert_eq!(geom(front), geom(back), "one shared plane chart");
}

fn generated_zero(brep: &BRep) -> FaceHandle {
    brep.faces
        .iter()
        .find(|(_, f)| {
            f.name
                .and_then(|n| brep.name_pool.get(n))
                .map(|n| n.role.category == RoleCategory::Generated && n.role.index == 0)
                .unwrap_or(false)
        })
        .map(|(h, _)| h)
        .expect("Generated/0 face")
}

#[test]
fn counter_clockwise_loop_lamina_is_flag_honest() {
    let mut sk = Sketch::new();
    rectangle(&mut sk, 0.0, 0.0, 1.0, 1.0, false);
    let (r, _) = build(sk);
    assert_eq!(r.brep.faces.len(), 2);
    let front = r.regions[0].face;
    assert_lamina_honest(&r.brep, front);
    assert_eq!(
        generated_zero(&r.brep),
        front,
        "Generated/0 is the front face"
    );
}

/// The tracer's direction (and hence which Euler-op face is the
/// enclosed one) follows the authored direction; the flags must not.
#[test]
fn clockwise_loop_lamina_is_flag_honest() {
    let mut sk = Sketch::new();
    rectangle(&mut sk, 0.0, 0.0, 1.0, 1.0, true);
    let (r, _) = build(sk);
    assert_eq!(r.brep.faces.len(), 2);
    let front = r.regions[0].face;
    assert_lamina_honest(&r.brep, front);
    assert_eq!(
        generated_zero(&r.brep),
        front,
        "Generated/0 is the front face"
    );
}

/// A region with a hole: the twin hole pcurves run the reversed span,
/// so the strict curve-vs-pcurve consistency check passes on both
/// faces.
#[test]
fn lamina_with_hole_is_flag_honest_and_pcurve_consistent() {
    let mut sk = Sketch::new();
    rectangle(&mut sk, 0.0, 0.0, 2.0, 2.0, false);
    rectangle(&mut sk, 0.5, 0.5, 1.5, 1.5, true);
    let (r, _) = build(sk);
    assert_eq!(r.brep.faces.len(), 2);
    assert_eq!(r.regions.len(), 1);
    assert_eq!(r.regions[0].holes.len(), 1);
    let front = r.regions[0].face;
    assert_lamina_honest(&r.brep, front);
    for (_, f) in &r.brep.faces {
        assert_eq!(f.inner_loops.len(), 1, "both faces carry the hole");
    }
}
