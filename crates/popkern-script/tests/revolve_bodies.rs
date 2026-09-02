//! PK-408: revolves recorded through the recorder serve VALID bodies —
//! a welded ring, a torus from a circle profile, a disc whose profile
//! touches the axis, a capped partial revolve — and a boolean on a
//! revolved ring serves. Before, every revolve replayed "served" and
//! failed Part validation (unwelded faces, inconsistent orientation).

use popkern_history::feature::AxisDesc;
use popkern_script::prelude::*;
use popkern_script::{Plane, RevolveEnd};

fn z_axis_at(x: f64) -> AxisDesc {
    AxisDesc::Literal {
        origin: [x, 0.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    }
}

fn counts(m: &mut Model) -> (usize, usize, usize) {
    let brep = m.to_brep().expect("brep");
    (brep.vertices.len(), brep.edges.len(), brep.faces.len())
}

#[test]
fn full_revolution_of_a_rectangle_is_a_valid_ring() {
    let mut m = Model::new();
    m.profile(Plane::XZ)
        .rect(0.2, 0.1)
        .at(0.6, 3.0)
        .sketch("ring_section")
        .expect("sketch")
        .revolve(&mut m, z_axis_at(0.0), RevolveEnd::FullCircle, "ring")
        .expect("revolve");
    assert_eq!(counts(&mut m), (4, 6, 4));
    let part = m.to_part().expect("ring validates");
    assert_eq!(part.face_count(), 4);
    // PK-409: the faces are exact (two cylinders, two planes) and the
    // mesh welds across every shared circle.
    let kinds = surface_kinds(part.brep());
    assert_eq!(
        kinds,
        ["Cylinder", "Cylinder", "Plane", "Plane"],
        "{kinds:?}"
    );
    m.assert_watertight().expect("watertight");
}

fn surface_kinds(brep: &popkern_core::topology::BRep) -> Vec<&'static str> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let mut kinds: Vec<&'static str> = brep
        .faces
        .values()
        .map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => "Plane",
            Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder",
            Some(SurfaceGeometry::Cone { .. }) => "Cone",
            Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
            Some(SurfaceGeometry::Torus { .. }) => "Torus",
            _ => "Nurbs",
        })
        .collect();
    kinds.sort();
    kinds
}

#[test]
fn circle_profile_revolves_to_a_torus_with_its_curvature() {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .circle(0.3)
        .at(1.0, 0.0)
        .sketch("tube_section")
        .expect("sketch")
        .revolve(
            &mut m,
            AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 1.0, 0.0],
            },
            RevolveEnd::FullCircle,
            "torus",
        )
        .expect("revolve");
    let brep = m.to_brep().expect("brep");
    // One periodic face, one seam circle + one profile circle, one vertex.
    assert_eq!(
        (brep.vertices.len(), brep.edges.len(), brep.faces.len()),
        (1, 2, 1)
    );
    let radii: Vec<f64> = brep
        .vertices
        .values()
        .map(|v| (v.position.x.powi(2) + v.position.z.powi(2)).sqrt())
        .collect();
    assert!(
        radii
            .iter()
            .all(|r| (r - 1.3).abs() < 1e-9 || (r - 0.7).abs() < 1e-9),
        "{radii:?}"
    );
    let part = m.to_part().expect("torus validates");
    assert_eq!(surface_kinds(part.brep()), ["Torus"]);
    m.assert_watertight().expect("watertight");
}

#[test]
fn profile_touching_the_axis_revolves_to_a_disc() {
    let mut m = Model::new();
    m.profile(Plane::XZ)
        .rect(0.2, 0.1)
        .at(0.1, 3.0)
        .sketch("disc_section")
        .expect("sketch")
        .revolve(&mut m, z_axis_at(0.0), RevolveEnd::FullCircle, "disc")
        .expect("revolve");
    assert_eq!(counts(&mut m), (2, 3, 3));
    let part = m.to_part().expect("disc validates");
    assert_eq!(surface_kinds(part.brep()), ["Cylinder", "Plane", "Plane"]);
    m.assert_watertight()
        .expect("watertight: the pole faces mesh closed");
}

#[test]
fn partial_revolution_is_capped_and_valid() {
    let mut m = Model::new();
    m.profile(Plane::XZ)
        .rect(0.2, 0.1)
        .at(0.6, 3.0)
        .sketch("arc_section")
        .expect("sketch")
        .revolve(
            &mut m,
            z_axis_at(0.0),
            RevolveEnd::Blind(90.0.into()),
            "quarter",
        )
        .expect("revolve");
    assert_eq!(counts(&mut m), (8, 12, 6));
    m.to_part().expect("quarter validates");
    m.assert_watertight()
        .expect("watertight: capped partial revolve");
}

#[test]
fn a_boolean_on_a_revolved_ring_serves() {
    let mut m = Model::new();
    let ring = m
        .profile(Plane::XZ)
        .rect(0.2, 0.1)
        .at(0.6, 3.0)
        .sketch("ring_section")
        .expect("sketch")
        .revolve(&mut m, z_axis_at(0.0), RevolveEnd::FullCircle, "ring")
        .expect("revolve");
    let cutter = m
        .tool(|m| {
            m.profile(Plane::xy_at(2.9))
                .circle(0.05)
                .at(0.6, 0.0)
                .sketch("cutter_sketch")?
                .extrude(m, End::Blind(0.3.into()))
        })
        .expect("cutter");
    ring.subtract(&mut m, cutter).expect("subtract serves");
    let part = m.to_part().expect("bored ring validates");
    assert!(
        part.face_count() > 4,
        "the bore adds faces: {}",
        part.face_count()
    );
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("revolve_bodies");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("bored_ring.pkn");
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("open");
    assert_eq!(again.replay().expect("replay").failures().count(), 0);
    assert_eq!(
        again.to_part().expect("part").face_count(),
        part.face_count()
    );
}

#[test]
fn slanted_line_revolves_to_an_exact_cone() {
    let mut m = Model::new();
    m.profile(Plane::XZ)
        .polygon(&[(0.0, 0.0), (0.5, 0.0), (0.0, 1.0)])
        .sketch("cone_section")
        .expect("sketch")
        .revolve(&mut m, z_axis_at(0.0), RevolveEnd::FullCircle, "cone")
        .expect("revolve");
    // Base disc + lateral cone with its apex; the on-axis span is dropped.
    assert_eq!(counts(&mut m), (2, 2, 2));
    let part = m.to_part().expect("cone validates");
    assert_eq!(surface_kinds(part.brep()), ["Cone", "Plane"]);
    m.assert_watertight()
        .expect("watertight: apex meshes closed");
}

/// STEP export of the exact ring passes the OCC oracle when the OCP
/// venv is present (skips otherwise).
#[test]
fn exact_ring_step_export_passes_the_occ_oracle() {
    let mut m = Model::new();
    m.profile(Plane::XZ)
        .rect(0.2, 0.1)
        .at(0.6, 3.0)
        .sketch("ring_section")
        .expect("sketch")
        .revolve(&mut m, z_axis_at(0.0), RevolveEnd::FullCircle, "ring")
        .expect("revolve");
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("revolve_bodies");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join("ring.step");
    m.assert_step_roundtrip(&step).expect("step round trip");
    let oracle = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let script = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tools/occ-oracle/occ_check.py");
    if oracle.exists() && script.exists() {
        let output = std::process::Command::new(&oracle)
            .arg(&script)
            .arg(&step)
            .output()
            .expect("oracle");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(report.contains("\"passed\": true"), "OCC oracle: {report}");
    } else {
        eprintln!("[pk409] OCC oracle absent; STEP round trip checked by popkern only");
    }
}
