//! Thermal contact interfaces (PK-1073): two blocks in series through
//! a contact conductance, against the exact series-resistance profile.

use super::*;
use crate::contact::ContactPatch;
use crate::mesh_test_helpers::box_mesh;
use crate::mesh_test_models::model_of_many;
use crate::model::SolidId;

fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-12,
        coarsest_dofs: 500,
        ..AnalysisOptions::default()
    }
}

/// Blocks [0, 20] and [20, 40] × [0, 10]² as two solids, faces 0..6 and
/// 6..12 (each box's x-min, x-max, then the rest).
fn two_blocks() -> AnalysisModel {
    model_of_many(
        vec![
            box_mesh([0.0; 3], [20.0, 10.0, 10.0]),
            box_mesh([20.0, 0.0, 0.0], [40.0, 10.0, 10.0]),
        ],
        6,
    )
}

/// The x = 20 square, wound with its normal +x (out of the first block).
fn interface() -> ContactPatch {
    let p = [
        [20.0, 0.0, 0.0],
        [20.0, 10.0, 0.0],
        [20.0, 10.0, 10.0],
        [20.0, 0.0, 10.0],
    ];
    ContactPatch {
        solids: (SolidId(0), SolidId(1)),
        triangles: vec![[p[0], p[1], p[2]], [p[0], p[2], p[3]]],
    }
}

fn problem(contacts: Vec<ThermalContact>) -> ConductionProblem {
    let mut faces = vec![FaceCondition::Insulated; 12];
    let fixed = |t0_c| FaceCondition::Fixed {
        t0_c,
        gradient: [0.0; 3],
    };
    faces[0] = fixed(100.0);
    faces[7] = fixed(20.0);
    ConductionProblem {
        conductivity: 167.0,
        faces,
        source_w: 0.0,
        body_sources: Vec::new(),
        hdt_c: None,
        contacts,
    }
}

#[test]
fn blocks_in_series_through_a_contact_have_the_series_resistance_profile() {
    let model = two_blocks();
    let hc = 5000.0;
    let contact = ThermalContact {
        patch: interface(),
        conductance: hc,
    };
    // Cells cutting the interface and cells whose planes hold it (the
    // patch on a grid plane is given to one cell; the other block's
    // fragment is across it).
    for h in [0.7, 0.5, 1.0, 2.0] {
        let r = solve_conduction(&model, &problem(vec![contact.clone()]), &options(h))
            .unwrap_or_else(|e| panic!("h {h}: {e}"));
        check(&model, &r, hc, h);
    }
}

fn check(model: &AnalysisModel, r: &ThermalResult, hc: f64, h: f64) {
    // Series resistances, K/W (k in W/(m·K), lengths in m).
    let (k, area, length) = (167.0, 1e-4, 0.02);
    let r_block = length / (k * area);
    let r_contact = 1.0 / (hc * area);
    let q = 80.0 / (2.0 * r_block + r_contact);
    let exact = |s: usize, x: f64| {
        let x_m = 1e-3 * x;
        if s == 0 {
            100.0 - q * x_m / (k * area)
        } else {
            20.0 + q * (0.04 - x_m) / (k * area)
        }
    };
    let mut worst = 0.0f64;
    for (s, temps) in r.vertex_temperature.iter().enumerate() {
        for (v, t) in model.solids[s].vertices.iter().zip(temps) {
            worst = worst.max((t - exact(s, v[0])).abs());
        }
    }
    assert!(worst <= 1e-6 * 80.0, "h {h}: worst vertex error {worst}");
    let into = -r.heat.face_out[0];
    let out = r.heat.face_out[7];
    assert!((into - q).abs() <= 1e-6 * q, "h {h}: into {into} vs {q}");
    assert!((out - q).abs() <= 1e-6 * q, "h {h}: out {out} vs {q}");
}

/// Control: without the contact the blocks pass no heat.
#[test]
fn without_the_contact_the_blocks_pass_no_heat() {
    let model = two_blocks();
    let r = solve_conduction(&model, &problem(Vec::new()), &options(0.7))
        .unwrap_or_else(|e| panic!("{e}"));
    assert!(r.heat.face_out[0].abs() < 1e-9, "{}", r.heat.face_out[0]);
}

/// A patch that lies inside one block and off the other declines: the
/// coupling would hold nothing.
#[test]
fn a_contact_patch_off_its_bodies_declines() {
    let model = two_blocks();
    let mut patch = interface();
    for t in &mut patch.triangles {
        for p in t {
            p[0] = 5.0;
        }
    }
    let contact = ThermalContact {
        patch,
        conductance: 5000.0,
    };
    let e = solve_conduction(&model, &problem(vec![contact]), &options(0.7)).err();
    assert!(
        matches!(e, Some(AnalysisError::Unsupported { .. })),
        "{e:?}"
    );
}
