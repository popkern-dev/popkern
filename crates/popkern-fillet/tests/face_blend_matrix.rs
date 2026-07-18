//! PK-142 stage 1 — face–face blend acceptance (SheetOnly).
//!
//! Fixture: a U-channel prism; its two inner walls face each other
//! across the slot and share no edge. The blend ball rolls in the
//! slot; the strip must come back as a valid open sheet whose
//! contact polylines lie on the two walls.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_fillet::{FaceBlendAttach, FaceBlendSpec, face_blend};

fn default_tol() -> Tolerance {
    Tolerance::default()
}

/// U-channel: outer 10×6, slot 4 wide × 4 deep from the top.
fn u_channel() -> BRep {
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(10.0, 6.0, 0.0),
        Point3::new(7.0, 6.0, 0.0),
        Point3::new(7.0, 2.0, 0.0),
        Point3::new(3.0, 2.0, 0.0),
        Point3::new(3.0, 6.0, 0.0),
        Point3::new(0.0, 6.0, 0.0),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), 5.0, default_tol()).unwrap()
}

/// The two inner slot walls: planes with normals ±X at x=3 and x=7.
fn slot_walls(brep: &BRep) -> (FaceHandle, FaceHandle) {
    let mut left = None;
    let mut right = None;
    for (fh, f) in brep.faces.iter() {
        let Some(sof) = f.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry else {
            continue;
        };
        if normal.x.abs() > 0.99 {
            if (origin.x - 3.0).abs() < 1e-9 {
                left = Some(fh);
            } else if (origin.x - 7.0).abs() < 1e-9 {
                right = Some(fh);
            }
        }
    }
    (
        left.expect("left slot wall"),
        right.expect("right slot wall"),
    )
}

#[test]
fn plane_plane_gap_sheet_blend() {
    let brep = u_channel();
    let (wall_a, wall_b) = slot_walls(&brep);
    // Slot width 4; parallel walls admit tangency at radius = gap/2.
    let result = face_blend(
        &brep,
        &FaceBlendSpec {
            face_a: wall_a,
            face_b: wall_b,
            radius: 2.0,
            attach: FaceBlendAttach::SheetOnly,
        },
        default_tol(),
    )
    .expect("plane-plane gap blend must complete");
    assert!(
        result.contacts_a.len() >= 2,
        "march must produce sections: {}",
        result.contacts_a.len()
    );
    // Contacts lie ON their walls.
    for p in &result.contacts_a {
        assert!((p.x - 3.0).abs() < 1e-6, "contact A off wall x=3: {:?}", p);
    }
    for p in &result.contacts_b {
        assert!((p.x - 7.0).abs() < 1e-6, "contact B off wall x=7: {:?}", p);
    }
    // The sheet is a structurally valid open sheet: one face, four
    // boundary edges, validate() clean under the open-sheet policy.
    assert_eq!(
        result.sheet.faces.len(),
        2,
        "two-sided strip sheet (front + back, the sheet_polygon model)"
    );
    assert_eq!(result.sheet.edges.len(), 4, "four boundary edges");
    result
        .sheet
        .validate()
        .expect("strip sheet must validate as an open sheet");
    // The strip surface midline stays medial: sample the surface
    // centre row and check equidistance from both walls to 1e-6.
    if let Some(sof) = result.sheet.faces.values().next().unwrap().surface.as_ref()
        && let SurfaceGeometry::Nurbs(n) = &sof.geometry
    {
        use popkern_core::geometry::SurfaceEval;
        let g = SurfaceGeometry::Nurbs(n.clone());
        let p = g.point_at(0.5, 0.5);
        assert!(
            (p.x - 5.0).abs() < 1e-6,
            "mid-strip point should sit on the slot midplane: {p:?}"
        );
    }
}

#[test]
fn gap_too_wide_declines_structured() {
    let brep = u_channel();
    let (wall_a, wall_b) = slot_walls(&brep);
    let err = face_blend(
        &brep,
        &FaceBlendSpec {
            face_a: wall_a,
            face_b: wall_b,
            radius: 1.0, // diameter 2 < gap 4 — parallel walls need gap/2
            attach: FaceBlendAttach::SheetOnly,
        },
        default_tol(),
    )
    .expect_err("gap 4 with ball diameter 2 must decline");
    let msg = format!("{err:?}");
    assert!(
        msg.contains("radius = gap/2") || msg.contains("exceeds the ball diameter"),
        "decline must name the gap class: {msg}"
    );
}

#[test]
fn touching_faces_decline_structured() {
    // Same face twice — zero gap; must route the caller to
    // fillet_edges instead of blending a face against itself.
    let brep = u_channel();
    let (wall_a, _) = slot_walls(&brep);
    let err = face_blend(
        &brep,
        &FaceBlendSpec {
            face_a: wall_a,
            face_b: wall_a,
            radius: 1.0,
            attach: FaceBlendAttach::SheetOnly,
        },
        default_tol(),
    )
    .expect_err("zero-gap pair must decline");
    let msg = format!("{err:?}");
    assert!(
        msg.contains("touch or intersect"),
        "decline must name the touching class: {msg}"
    );
}
