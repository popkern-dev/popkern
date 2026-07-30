//! PK-142 — cylinder-plane face–face blend: a vertical rod standing
//! near a tilted wall. Covers the two remaining acceptance classes:
//! curved-support contact curves (G1 sampled against both supports)
//! and a gap that varies along the march (the ball loses two-contact
//! where the gap exceeds the diameter, terminating the strip inside
//! the face trims).

use popkern_bool::union;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_fillet::{FaceBlendAttach, FaceBlendSpec, face_blend};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Base plate + leaning wall + vertical rod. The wall's east face
/// tilts away from the rod as z rises, so the wall–rod gap grows both
/// with height and around the rod's circumference.
fn rod_near_wall() -> BRep {
    let base = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(12.0, 0.0, 0.0),
            Point3::new(12.0, 10.0, 0.0),
            Point3::new(0.0, 10.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        2.0,
        tol(),
    )
    .expect("base builds");
    // Wall profile in the xz plane at y = 0, extruded along +y: east
    // face leans from x = 2.5 (z = 1) to x = 1.7 (z = 14).
    let wall = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(2.5, 0.0, 1.0),
            Point3::new(1.7, 0.0, 14.0),
            Point3::new(0.0, 0.0, 14.0),
        ],
        Vec3::new(0.0, 1.0, 0.0),
        10.0,
        tol(),
    )
    .expect("wall builds");
    let rod = {
        let cyl = BRep::cylinder(1.5, 12.0, tol()).expect("rod builds");
        let t = popkern_core::geometry::transform::Transform3::translation(5.3, 5.0, 1.0);
        cyl.transform(&t).expect("rod transform")
    };
    let body = union(base, wall, tol(), None).expect("base+wall union");
    union(body, rod, tol(), None).expect("body+rod union")
}

/// The wall's leaning east face and the rod's cylinder surface.
fn blend_pair(brep: &BRep) -> (FaceHandle, FaceHandle) {
    let mut wall = None;
    let mut rod = None;
    for (fh, f) in brep.faces.iter() {
        let Some(sof) = f.surface.as_ref() else {
            continue;
        };
        match &sof.geometry {
            SurfaceGeometry::Plane { normal, .. } => {
                if normal.x > 0.9 && normal.z.abs() > 0.01 {
                    wall = Some(fh);
                }
            }
            SurfaceGeometry::Cylinder { radius, .. } if (radius - 1.5).abs() < 1e-9 => {
                rod = Some(fh);
            }
            _ => {}
        }
    }
    (wall.expect("leaning wall face"), rod.expect("rod face"))
}

#[test]
fn cylinder_plane_gap_sheet_blend_g1() {
    let body = rod_near_wall();
    let (wall, rod) = blend_pair(&body);
    let radius = 0.75;
    let result = face_blend(
        &body,
        &FaceBlendSpec {
            face_a: wall,
            face_b: rod,
            radius,
            // Constant law: PK-162's linear radius law is exercised by
            // the fillet crate's own variable-radius gates.
            radius_end: None,
            attach: FaceBlendAttach::SheetOnly,
        },
        tol(),
    )
    .expect("cylinder-plane sheet blend must complete");
    assert!(
        result.contacts_a.len() >= 3,
        "march must produce a real span: {} sections",
        result.contacts_a.len()
    );
    result
        .sheet
        .validate()
        .expect("cylinder-plane strip sheet must validate");

    // G1 against both supports, sampled at every section: the ball
    // centre sits along each support normal at the contact.
    let axis_base = Point3::new(5.3, 5.0, 1.0);
    let g1_bar = 5.0e-2;
    for (pa, pb) in result.contacts_a.iter().zip(result.contacts_b.iter()) {
        let centre_dir_a = Vec3::new(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z);
        // Wall normal (unit): the leaning plane's outward normal.
        let n_wall = Vec3::new(13.0, 0.0, 0.8).normalized().expect("wall normal");
        let along = centre_dir_a.normalized().expect("section direction");
        // contact-to-contact direction need not be the normal (the
        // section is an arc), so measure each side independently:
        // plane side — the contact must be the foot of the centre.
        let _ = along;
        // Reconstruct the centre from the plane contact + normal.
        let centre = Point3::new(
            pa.x + n_wall.x * radius,
            pa.y + n_wall.y * radius,
            pa.z + n_wall.z * radius,
        );
        // Cylinder side: the ball sits in the gap OUTSIDE the rod, so
        // its centre lies radially outward of the rod contact (axis
        // distance R + r).
        let foot = Point3::new(axis_base.x, axis_base.y, pb.z);
        let radial = Vec3::new(pb.x - foot.x, pb.y - foot.y, 0.0)
            .normalized()
            .expect("radial direction");
        let expect_b = Point3::new(
            pb.x + radial.x * radius,
            pb.y + radial.y * radius,
            pb.z + radial.z * radius,
        );
        let mismatch = centre.distance_to(&expect_b);
        assert!(
            mismatch < g1_bar * radius,
            "G1 violation: plane-implied centre {centre:?} vs cylinder-implied \
             {expect_b:?} (mismatch {mismatch:.4})"
        );
    }

    // Variable gap: the strip must TERMINATE where the ball loses
    // two-contact — well inside the wall's y-extent [0, 10].
    let (mut y_min, mut y_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in &result.contacts_a {
        y_min = y_min.min(p.y);
        y_max = y_max.max(p.y);
    }
    assert!(
        y_min > 0.5 && y_max < 9.5,
        "strip must end where the gap exceeds the ball diameter, inside the \
         trim: y span [{y_min:.2}, {y_max:.2}]"
    );
}
