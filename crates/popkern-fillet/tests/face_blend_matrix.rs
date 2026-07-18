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

/// Truncated-V notch: outer 10×6, tilted walls converging to a
/// width-1 floor at y=2 from a width-3 opening at the top. The ball
/// position is fully determined by the wall angle (no line-contact
/// ambiguity), so this fixture pins TrimBoth numerically.
fn v_notch() -> BRep {
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(10.0, 6.0, 0.0),
        Point3::new(6.5, 6.0, 0.0),
        Point3::new(5.5, 2.0, 0.0),
        Point3::new(4.5, 2.0, 0.0),
        Point3::new(3.5, 6.0, 0.0),
        Point3::new(0.0, 6.0, 0.0),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), 5.0, default_tol()).unwrap()
}

/// The two tilted notch walls, identified by their outward normals
/// (±x dominant, tilted in y).
fn notch_walls(brep: &BRep) -> (FaceHandle, FaceHandle) {
    let mut left = None;
    let mut right = None;
    for (fh, f) in brep.faces.iter() {
        let Some(sof) = f.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Plane { normal, .. } = &sof.geometry else {
            continue;
        };
        // The notch walls are the only faces tilted in y with a
        // dominant x normal (outer ±x faces have normal.y == 0).
        if normal.z.abs() > 0.01 || normal.x.abs() < 0.5 || normal.y.abs() < 0.2 {
            continue;
        }
        if normal.x < 0.0 {
            right = Some(fh); // wall from (6.5,6) to (5.5,2): outward -x
        } else {
            left = Some(fh); // wall from (4.5,2) to (3.5,6): outward +x
        }
    }
    (
        left.expect("left notch wall"),
        right.expect("right notch wall"),
    )
}

#[test]
fn trim_both_through_notch_declines_open_ends() {
    // The V-notch is a THROUGH feature: the gap pocket is open at
    // both prism ends, so the strip end arcs land over void — the
    // closed-pocket trim surgery must decline by name (the
    // material-adding attach is the recorded follow-up). The blend
    // GEOMETRY (sheet mode) still completes and is pinned
    // analytically below.
    let brep = v_notch();
    let (wall_a, wall_b) = notch_walls(&brep);
    let err = face_blend(
        &brep,
        &FaceBlendSpec {
            face_a: wall_a,
            face_b: wall_b,
            radius: 1.0,
            attach: FaceBlendAttach::TrimBoth,
        },
        default_tol(),
    )
    .expect_err("through-notch TrimBoth must decline");
    let msg = format!("{err}");
    assert!(
        msg.contains("end-arc support"),
        "decline must name the end-arc mechanism: {msg}"
    );
}

#[test]
fn sheet_only_v_notch_contacts_analytic() {
    // Contact geometry is determined by the wall tilt alone: the
    // r = 1 ball centre sits at y = sqrt(17), contacts at
    // y = sqrt(17) - 1/sqrt(17) ≈ 3.88 — pinned analytically.
    let brep = v_notch();
    let (wall_a, wall_b) = notch_walls(&brep);
    let result = face_blend(
        &brep,
        &FaceBlendSpec {
            face_a: wall_a,
            face_b: wall_b,
            radius: 1.0,
            attach: FaceBlendAttach::SheetOnly,
        },
        default_tol(),
    )
    .expect("V-notch sheet blend must complete");
    assert!(result.contacts_a.len() >= 2);
    let y_expect = 17.0_f64.sqrt() - 1.0 / 17.0_f64.sqrt();
    for p in result.contacts_a.iter().chain(result.contacts_b.iter()) {
        assert!(
            (p.y - y_expect).abs() < 1e-6,
            "contact height off the analytic tangency: {p:?} vs y={y_expect}"
        );
    }
    result
        .sheet
        .validate()
        .expect("V-notch strip sheet must validate");
}

/// Deep channel for the parallel TrimBoth leg: slot 4 wide × 6 deep
/// in a 10×7 plate extruded 8. The canonical mid-extent line contact
/// (y = 4) leaves the r = 2 strip clear of the floor at y = 1, and
/// the z extent (8) dominates so the march runs along the extrusion.
fn deep_channel() -> BRep {
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(10.0, 7.0, 0.0),
        Point3::new(7.0, 7.0, 0.0),
        Point3::new(7.0, 1.0, 0.0),
        Point3::new(3.0, 1.0, 0.0),
        Point3::new(3.0, 7.0, 0.0),
        Point3::new(0.0, 7.0, 0.0),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), 8.0, default_tol()).unwrap()
}

#[test]
fn trim_both_deep_channel_declines_open_ends() {
    // The deep channel is open at both extrusion ends — same
    // open-ended class as the through notch, via the parallel
    // (line-contact) march. The canonical placement is still pinned:
    // line contacts at the wall mid-extent height, half-circle
    // sections bowing toward the floor (the closed side).
    let brep = deep_channel();
    let (wall_a, wall_b) = slot_walls(&brep);
    let err = face_blend(
        &brep,
        &FaceBlendSpec {
            face_a: wall_a,
            face_b: wall_b,
            radius: 2.0,
            attach: FaceBlendAttach::TrimBoth,
        },
        default_tol(),
    )
    .expect_err("open-ended channel TrimBoth must decline");
    let msg = format!("{err}");
    assert!(
        msg.contains("end-arc support"),
        "decline must name the end-arc mechanism: {msg}"
    );

    // Sheet mode: canonical mid-extent line contact + floor-bowed
    // half-circle sections.
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
    .expect("deep-channel sheet blend must complete");
    for p in result.contacts_a.iter().chain(result.contacts_b.iter()) {
        assert!(
            (p.y - 4.0).abs() < 1e-6,
            "parallel line contact must sit at the mid-extent height: {p:?}"
        );
    }
    // The section is a genuine half-cylinder bowing toward the floor
    // (y-), not a flat chord diaphragm: the sheet surface's mid-v row
    // must dip below the contact height by the radius.
    let sof = result
        .sheet
        .faces
        .values()
        .next()
        .and_then(|f| f.surface.as_ref())
        .expect("sheet has a surface");
    if let SurfaceGeometry::Nurbs(n) = &sof.geometry {
        use popkern_core::geometry::SurfaceEval;
        let g = SurfaceGeometry::Nurbs(n.clone());
        let p = g.point_at(0.5, 0.5);
        assert!(
            (p.y - 2.0).abs() < 1e-6,
            "half-circle section must bow to the floor side: {p:?}"
        );
        assert!((p.x - 5.0).abs() < 1e-6, "bow stays on the midplane: {p:?}");
    }
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
