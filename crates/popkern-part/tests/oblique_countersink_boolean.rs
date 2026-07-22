//! PK-166: oblique countersink boolean — a coaxial drafted-frustum
//! subtract whose cone crosses a TILTED rest plane (the laptop-stand
//! round-seat mouth class). The perpendicular-mouth class is gated in
//! cone_band_dressup_tests; this fixture pins the oblique class:
//! cone×oblique-plane = an exact ellipse (PK-154), and the boolean
//! must serve the cone lead-in band, the bore wall trimmed at the
//! crossing circle, and the rest-plane hole bounded by the cone-rim
//! ellipse. Measured at authoring: the cone face was never built —
//! a flat annular band shipped instead and OCC flagged
//! BRepCheck_IntersectingWires.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_part::{Dir, FaceQuery, Plane, Profile, SurfaceKind};

/// Vertical-ray parity containment against the default-quality mesh
/// (the cone_band_dressup instrument, ray nudged off seam columns).
fn contains(part: &popkern_part::Part, q: (f64, f64, f64)) -> bool {
    let m = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let q = (q.0 + 1.03e-7, q.1 + 1.71e-7, q.2);
    let mut crossings = 0;
    for t in &m.triangles {
        let (a, b, c) = (
            m.vertices[t[0]].position,
            m.vertices[t[1]].position,
            m.vertices[t[2]].position,
        );
        let (e1x, e1y) = (b.x - a.x, b.y - a.y);
        let (e2x, e2y) = (c.x - a.x, c.y - a.y);
        let det = e1x * e2y - e1y * e2x;
        if det.abs() < 1e-14 {
            continue;
        }
        let (rx, ry) = (q.0 - a.x, q.1 - a.y);
        let u = (rx * e2y - ry * e2x) / det;
        let v = (e1x * ry - e1y * rx) / det;
        if u < -1e-12 || v < -1e-12 || u + v > 1.0 + 1e-12 {
            continue;
        }
        let z = a.z + u * (b.z - a.z) + v * (c.z - a.z);
        if z > q.2 {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}

/// The fixture: a 60×60×40 plate whose top is shaved to a 12° oblique
/// plane (through (0,0,30/cos12°), normal (0,−sin12°,cos12°)), a
/// vertical through-bore r=6, then a coaxial countersink frustum
/// (r 5.4 at z 24, flaring 8.5° upward) whose cone crosses the bore
/// wall in the exact circle at z≈28.0 and the tilted face in an
/// ellipse strictly outside the bore mouth.
fn oblique_countersunk_plate() -> popkern_part::Part {
    let plate = Profile::on(Plane::XY)
        .rect(60.0, 60.0)
        .extrude(40.0)
        .expect("plate");
    let shave = Profile::on(Plane::xy_at(30.0))
        .rect(200.0, 200.0)
        .extrude(100.0)
        .expect("shave slab")
        .rotate_x(12.0_f64.to_radians())
        .expect("tilt");
    let plate = plate.subtract(&shave).expect("tilted top");
    let bore = Profile::on(Plane::xy_at(-1.0))
        .circle(6.0)
        .extrude(45.0)
        .expect("bore");
    let plate = plate.subtract(&bore).expect("bore cut");
    let sink = Profile::on(Plane::xy_at(24.0))
        .circle(5.4)
        .extrude(20.0)
        .expect("sink cyl")
        .draft(-8.5_f64.to_radians(), Dir::Z, |f: FaceQuery| {
            f.kind(SurfaceKind::Cylinder)
        })
        .expect("sink draft");
    plate.subtract(&sink).expect("oblique countersink")
}

#[test]
fn oblique_countersink_serves_cone_band_and_watertight_shell() {
    let cut = oblique_countersunk_plate();
    assert_eq!(
        cut.open_edge_count().expect("mesh"),
        0,
        "oblique countersink opens the shell"
    );
    // The lead-in must be an analytic CONE face — the measured
    // authoring defect shipped a flat annulus instead.
    let cone_faces = cut
        .brep()
        .faces
        .values()
        .filter(|f| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cone { .. }))
        })
        .count();
    assert_eq!(
        cone_faces, 1,
        "oblique countersink must serve exactly one cone lead-in face"
    );
    // 4 plate sides + bottom (with bore hole) + tilted top (with
    // cone-rim ellipse hole) + bore wall + cone band.
    assert_eq!(cut.face_count(), 8, "oblique countersink face census");

    // Material probes. Cone radius r(z) = 5.4 + (z−24)·tan 8.5°;
    // crossing circle with the r=6 bore at z ≈ 28.01.
    // Above the crossing the cone has eaten past the bore wall:
    assert!(
        !contains(&cut, (6.15, 0.0, 29.5)),
        "cone void missing above the crossing circle"
    );
    // Below the crossing the bore wall is intact (cone is narrower):
    assert!(
        contains(&cut, (6.3, 0.0, 26.0)),
        "bore wall breached below the crossing circle"
    );
    // Bore void itself:
    assert!(!contains(&cut, (5.0, 0.0, 20.0)), "bore void missing");
    // Outside the cone mouth on the tilted face, material remains
    // (tilted top at x-axis: z = 30/cos12° ≈ 30.68):
    assert!(
        contains(&cut, (8.0, 0.0, 30.0)),
        "material outside the cone mouth vanished"
    );
    cut.brep()
        .validate()
        .expect("oblique countersink validates");
}

#[test]
fn oblique_countersink_exports_clean_step() {
    let cut = oblique_countersunk_plate();
    let dir = std::path::Path::new("../../examples/output/fixtures");
    std::fs::create_dir_all(dir).expect("fixture output dir");
    let path = dir.join("oblique_countersink.step");
    cut.export_step(&path).expect("oblique countersink exports");
    cut.assert_step_roundtrip(&path)
        .expect("oblique countersink round-trips through STEP");
}

/// The stand-parameter variant: 18° rest tilt, BLIND bore (seat
/// floor at z 70), Ø70.4 seat, sink r 34.6 at z 82 flaring 8.5° so
/// the crossing circle sits at z 86 below the whole mouth — the
/// exact laptop-stand round-seat class (PK-166). The blind floor
/// adds a third ring (the floor rim) to the nested-hole assignment.
fn stand_class_countersunk_plate() -> popkern_part::Part {
    let plate = Profile::on(Plane::XY)
        .rect(200.0, 200.0)
        .extrude(140.0)
        .expect("plate");
    // Rotate about the origin FIRST, then lift — the rest plane is
    // z = 100 + y·tan(18°) exactly (the laptop-stand construction).
    let shave = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("shave slab")
        .rotate_x(18.0_f64.to_radians())
        .expect("tilt")
        .translate(0.0, 0.0, 100.0)
        .expect("lift");
    let plate = plate.subtract(&shave).expect("tilted top");
    let seat = Profile::on(Plane::xy_at(70.0))
        .circle(35.2)
        .extrude(60.0)
        .expect("seat");
    let plate = plate.subtract(&seat).expect("seat cut");
    let sink = Profile::on(Plane::xy_at(82.0))
        .circle(34.6)
        .extrude(36.0)
        .expect("sink cyl")
        .draft(-8.5_f64.to_radians(), Dir::Z, |f: FaceQuery| {
            f.kind(SurfaceKind::Cylinder)
        })
        .expect("sink draft");
    plate.subtract(&sink).expect("stand-class countersink")
}

#[test]
#[ignore = "PK-166 pin: the stand-class sink subtract's cellular walk breaks on the bore \
cylinder's chart (degree imbalance on its two closed periodic records — the crossing \
circle and the mouth ellipse — at the periodic seam), the cone fragment's boundary \
conflicts mint duplicate edges instead of sharing, and a phantom outer-1/inner-1 \
fragment ships on the tilted face: 741 open mesh edges, no cone face. Un-ignore when \
the cellular chart handles multi-closed-record cylinder faces."]
fn stand_class_oblique_countersink_serves_cone_band() {
    let cut = stand_class_countersunk_plate();
    assert_eq!(
        cut.open_edge_count().expect("mesh"),
        0,
        "stand-class countersink opens the shell"
    );
    let cone_faces = cut
        .brep()
        .faces
        .values()
        .filter(|f| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cone { .. }))
        })
        .count();
    assert_eq!(
        cone_faces, 1,
        "stand-class countersink must serve exactly one cone lead-in face"
    );
    // Probes: crossing circle at z 86 (r_cone(86) = 34.6 + 4·tan8.5°
    // = 35.2). Above it the cone has eaten past the bore wall:
    assert!(
        !contains(&cut, (35.5, 0.0, 88.5)),
        "cone void missing above the crossing circle"
    );
    // Below the crossing the bore wall is intact:
    assert!(
        contains(&cut, (35.5, 0.0, 84.0)),
        "bore wall breached below the crossing circle"
    );
    // Seat void and floor:
    assert!(!contains(&cut, (0.0, 0.0, 80.0)), "seat void missing");
    assert!(contains(&cut, (0.0, 0.0, 60.0)), "seat floor missing");
    cut.brep()
        .validate()
        .expect("stand-class countersink validates");
}
