//! The adapter on real kernel BReps (PK-948). These are the kernel
//! suite's smoke set for the BRep -> `AnalysisModel` seam: a kernel
//! change to tessellation, measurement or primitive orientation that
//! breaks what the solver reads fails here, without waiting for the
//! physics suite.

use popkern_analysis::{AnalysisModel, FaceId, Provenance, Surface};
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::TessellationOptions;

use super::*;

fn provenance() -> Provenance {
    Provenance::Popkern {
        document: "adapter_tests".into(),
        revision: 0,
    }
}

fn model(breps: &[BRep]) -> AnalysisModel {
    model_from_breps(
        breps,
        &TessellationOptions::display(),
        &FeatureTolerances::default(),
        provenance(),
    )
    .unwrap_or_else(|e| panic!("{e}"))
}

/// Volume enclosed by a solid's welded mesh (positive when outward).
fn mesh_volume(model: &AnalysisModel, solid: usize) -> f64 {
    let s = &model.solids[solid];
    s.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| s.vertices[i as usize]);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

#[test]
fn a_unit_cube_has_six_outward_unit_planes() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let m = model(&[cube]);
    assert_eq!(m.solids.len(), 1);
    assert_eq!(m.faces.len(), 6);
    let exact = m.solids[0].exact_volume.expect("exact volume");
    assert!((exact - 1.0).abs() < 1e-12, "{exact}");
    assert!((mesh_volume(&m, 0) - 1.0).abs() < 1e-12);
    // The cube spans [0, 1] or [-0.5, 0.5]: read its centre off the mesh.
    let n = m.solids[0].vertices.len() as f64;
    let centre = m.solids[0].vertices.iter().fold([0.0; 3], |c, v| {
        [c[0] + v[0] / n, c[1] + v[1] / n, c[2] + v[2] / n]
    });
    for face in &m.faces {
        assert!((face.area - 1.0).abs() < 1e-12, "{face:?}");
        let Surface::Plane { normal, offset } = face.surface else {
            panic!("{face:?} is not a plane");
        };
        let at_centre = normal[0] * centre[0] + normal[1] * centre[1] + normal[2] * centre[2];
        assert!(
            offset - at_centre > 0.4,
            "face {:?}: normal {normal:?} does not point outward",
            face.id
        );
    }
}

#[test]
fn curved_primitives_mesh_close_to_their_exact_volumes() {
    let tol = Tolerance::default();
    // Turned off the world axes, so a chart or seam assumption tied to
    // +Z shows here.
    let turn = Transform3::rotation_x(0.7).compose(&Transform3::rotation_z(1.9));
    let cylinder = BRep::cylinder(0.5, 2.0, tol)
        .and_then(|c| c.transform(&turn))
        .expect("cylinder");
    let sphere = BRep::sphere(0.5, tol)
        .and_then(|s| s.transform(&turn))
        .expect("sphere");
    for (name, brep, kinds) in [
        ("cylinder", cylinder, ["Cylinder", "Plane"]),
        ("sphere", sphere, ["Sphere", "Sphere"]),
    ] {
        let m = model(&[brep]);
        let exact = m.solids[0].exact_volume.expect("exact volume");
        let meshed = mesh_volume(&m, 0);
        assert!(
            (meshed - exact).abs() < 0.01 * exact,
            "{name}: mesh volume {meshed} vs exact {exact}"
        );
        for face in &m.faces {
            let kind = format!("{:?}", face.surface);
            assert!(
                kinds.iter().any(|k| kind.starts_with(k)),
                "{name}: unexpected surface {kind}"
            );
            assert!(face.folded_area.unwrap_or(0.0) < 1e-9, "{name}: {face:?}");
        }
    }
}

#[test]
fn a_bored_block_is_watertight_and_keeps_its_volume() {
    let tol = Tolerance::default();
    let block = BRep::unit_cube(tol).expect("cube");
    let n = block.vertices.len() as f64;
    let c = block.vertices.values().fold([0.0; 3], |c, v| {
        [
            c[0] + v.position.x / n,
            c[1] + v.position.y / n,
            c[2] + v.position.z / n,
        ]
    });
    let bore = BRep::cylinder(0.2, 2.0, tol)
        .and_then(|b| b.transform(&Transform3::translation(c[0], c[1], c[2] - 1.0)))
        .expect("bore");
    let bored = boolean(BoolOperation::Subtraction, block, bore, tol, None).expect("bore");
    let m = model(&[bored]);
    let want = 1.0 - std::f64::consts::PI * 0.04;
    let exact = m.solids[0].exact_volume.expect("exact volume");
    assert!((exact - want).abs() < 1e-9, "exact {exact} vs {want}");
    let meshed = mesh_volume(&m, 0);
    assert!((meshed - want).abs() < 0.005 * want, "mesh {meshed}");
}

#[test]
fn named_faces_carry_the_namer_and_map_back_to_their_source() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let cylinder = BRep::cylinder(0.25, 1.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(3.0, 0.0, 0.0)))
        .expect("cylinder");
    let breps = [cube, cylinder];
    let named = model_from_breps_named(
        &breps,
        &TessellationOptions::display(),
        &FeatureTolerances::default(),
        provenance(),
        &|brep, face| (brep == 0).then(|| format!("cube.{face:?}")),
    )
    .expect("model");
    let m = &named.model;
    assert_eq!(m.solids.len(), 2);
    assert_eq!(named.face_sources.len(), m.faces.len());
    for face in &m.faces {
        let (brep, handle) = named.face_sources[face.id.0 as usize];
        assert!(breps[brep].faces.get(handle).is_some());
        match brep {
            0 => assert_eq!(face.name, Some(format!("cube.{handle:?}"))),
            _ => assert_eq!(face.name, None),
        }
    }
    // The cube's six faces come first, then the cylinder's.
    assert!(m.faces[..6].iter().all(|f| f.solid == m.solids[0].id));
    assert_eq!(m.faces[6].id, FaceId(6));
    // Naming changes nothing else.
    let mut unnamed = named.model.clone();
    for face in &mut unnamed.faces {
        face.name = None;
    }
    assert_eq!(model(&breps), unnamed);
}
