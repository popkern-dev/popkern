//! Operand material senses read off point classifications beside the
//! faces, with the volume sign as the fallback.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::MaterialSense;

use super::{material_senses, probe_votes};

/// A 10 × 10 prism `thickness` tall. `BRep::prism` folds its normals
/// INTO the material.
fn slab(thickness: f64, tol: Tolerance) -> BRep {
    BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
            Point3::new(0.0, 10.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        thickness,
        tol,
    )
    .expect("slab")
}

/// Whether `f` is one of the slab's caps (a plane facing ±z).
fn is_cap(b: &BRep, f: FaceHandle) -> bool {
    b.faces[f].surface.as_ref().is_some_and(
        |s| matches!(s.geometry, SurfaceGeometry::Plane { normal, .. } if normal.z.abs() > 0.5),
    )
}

/// A plate 1e-4 thick, under the probe offset (1e-5 of its 14 box
/// diagonal), probed caps first: beside a cap on either side lies
/// outside, so the caps cast no vote and the walls decide, as the
/// volume sign of the same prism built 1 thick says. A one-sided probe
/// stepped through both inward-folded caps and voted them outward.
#[test]
fn a_plate_thinner_than_the_probe_offset_keeps_its_sense() {
    let tol = Tolerance::default();
    let thick = slab(1.0, tol);
    let seed = thick.faces.keys().min().expect("faces");
    let want = thick.component_material_sense(seed) == Some(MaterialSense::FoldedOutward);
    assert!(!want, "prisms fold inward");
    let plate = slab(1.0e-4, tol);
    let boxes = crate::intersect::compute_face_aabbs(&plate);
    let (caps, walls): (Vec<FaceHandle>, Vec<FaceHandle>) =
        plate.faces.keys().partition(|&f| is_cap(&plate, f));
    assert_eq!((caps.len(), walls.len()), (2, 4));
    let caps_first: Vec<FaceHandle> = caps.into_iter().chain(walls).collect();
    assert_eq!(probe_votes(&plate, &caps_first, &boxes, &tol), [want, want]);
    let senses = material_senses(&plate, 0, &boxes, &tol).expect("decided");
    assert!(
        senses.values().all(|&outward| outward == want),
        "{senses:?}"
    );
}
