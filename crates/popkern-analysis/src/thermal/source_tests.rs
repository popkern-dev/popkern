//! Heat sources on named bodies (Stage 1b, PH-049).

use super::from_setup::solve_thermal;
use super::*;
use crate::face_ref::Tolerance;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::model::{FaceInfo, SolidId};
use crate::thermal_setup::{
    Convection, ConvectionFaces, HeatTransfer, ThermalSetup, ThermalSetupError, VolumeSource,
};

/// Two separate boxes, `heater` (faces 0..5) and `spacer` (6..11).
fn two_bodies() -> AnalysisModel {
    let mut m = model_of(box_mesh([0.0; 3], [6.0, 4.0, 4.0]), 6);
    m.solids[0].name = Some("heater".to_string());
    let mut second = box_mesh([10.0, 0.0, 0.0], [16.0, 4.0, 4.0]);
    second.id = SolidId(1);
    second.name = Some("spacer".to_string());
    for f in &mut second.triangle_face {
        f.0 += 6;
    }
    for f in 0..6u32 {
        let mut info: FaceInfo = m.faces[f as usize].clone();
        info.id = FaceId(6 + f);
        info.solid = SolidId(1);
        m.faces.push(info);
    }
    m.solids.push(second);
    m
}

fn setup(body: &str) -> ThermalSetup {
    ThermalSetup {
        material: "Aluminium 6061-T6".to_string(),
        tolerance: Tolerance::default(),
        faces: Vec::new(),
        fixed: Vec::new(),
        flux: Vec::new(),
        convection: vec![Convection {
            faces: ConvectionFaces::AllOther,
            h: HeatTransfer::Value(10.0),
            ambient_c: 20.0,
            emittance: None,
        }],
        sources: vec![VolumeSource {
            power_w: 2.0,
            body: Some(body.to_string()),
        }],
        transient: None,
        options: AnalysisOptions {
            target_cell_size: Some(0.8),
            cg_rel_tol: 1e-12,
            ..AnalysisOptions::default()
        },
    }
}

#[test]
fn a_source_in_one_body_leaves_through_that_body_only() {
    let model = two_bodies();
    let r = solve_thermal(&model, &setup("heater"), &[]).unwrap_or_else(|e| panic!("{e}"));
    let out = |faces: std::ops::Range<usize>| -> f64 { faces.map(|f| r.heat.face_out[f]).sum() };
    assert!(
        (out(0..6) - 2.0).abs() <= 1e-9 * 2.0,
        "heater {}",
        out(0..6)
    );
    assert!(out(6..12).abs() <= 1e-12, "spacer {}", out(6..12));
    assert_eq!(r.heat.source, 2.0);
    // The spacer stays at ambient; the heater is warmer.
    let spacer = &r.vertex_temperature[1];
    assert!(
        spacer.iter().all(|t| (t - 20.0).abs() <= 1e-9),
        "{spacer:?}"
    );
    assert!(r.vertex_temperature[0].iter().all(|t| *t > 20.0));
}

#[test]
fn an_unknown_body_declines_by_name() {
    let model = two_bodies();
    match solve_thermal(&model, &setup("boiler"), &[]) {
        Err(AnalysisError::InvalidThermalSetup(ThermalSetupError::UnknownBody { name, known })) => {
            assert_eq!(name, "boiler");
            assert_eq!(known, ["heater", "spacer"]);
        }
        other => panic!("{other:?}"),
    }
}
