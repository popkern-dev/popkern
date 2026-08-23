use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// A unit-square sheet with authored geometry (plane + line curves),
/// the state the history sketch path hands to extrude.
fn square_sheet() -> (BRep, crate::topology::FaceHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
    brep.attach_geometry().unwrap();
    brep.attach_planar_pcurves().unwrap();
    (brep, mef.face)
}

/// PK-316: extrude authors its own geometry — the result validates
/// WITHOUT any post-pass, every face carries a distinct correct
/// surface, and every edge carries a curve. (Before: every mef face
/// inherited the base plane and new edges had no curves — the history
/// replay path shipped validate-red bodies.)
#[test]
fn extrude_authors_geometry_and_validates() {
    let (mut brep, face) = square_sheet();
    let result = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();

    brep.validate()
        .expect("extruded body validates with no post-pass");

    for (eh, ed) in &brep.edges {
        assert!(ed.curve.is_some(), "edge {eh:?} must carry a curve");
    }
    // Side faces: four distinct planes, none of them z-normal.
    assert_eq!(result.side_faces.len(), 4);
    for &f in &result.side_faces {
        let sof = brep.faces[f].surface.as_ref().expect("side surface");
        let SurfaceGeometry::Plane { normal, .. } = &sof.geometry else {
            panic!("straight-edge side face must be planar");
        };
        assert!(
            normal.z.abs() < 1e-9,
            "side face normal {normal:?} must be lateral, not the base plane's"
        );
    }
    // Cap: the base plane translated to z = 1.
    let cap = brep.faces[result.cap_face].surface.as_ref().unwrap();
    let SurfaceGeometry::Plane { origin, normal, .. } = &cap.geometry else {
        panic!("cap must be planar");
    };
    assert!(
        (origin.z - 1.0).abs() < 1e-9,
        "cap origin at z=1, got {origin:?}"
    );
    assert!(normal.z.abs() > 0.99, "cap normal along z, got {normal:?}");
}

/// The full-validation screen is clean too (self-intersection etc.).
#[test]
fn extruded_sheet_passes_full_validation() {
    let (mut brep, face) = square_sheet();
    brep.extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();
    let findings = brep.validate_all();
    assert!(
        findings.is_empty(),
        "full validation clean, got {findings:?}"
    );
}
