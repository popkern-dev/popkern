use super::*;
use crate::mesh_test_helpers::{box_mesh, model_of, tube};
use crate::thermal::solve_conduction;

fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-12,
        coarsest_dofs: 500,
        ..AnalysisOptions::default()
    }
}

#[test]
fn a_linear_field_is_reproduced_through_the_caps() {
    let model = model_of(box_mesh([0.0; 3], [6.0, 4.0, 3.0]), 6);
    let g = [2.0, -1.0, 0.5];
    let problem = ConductionProblem {
        conductivity: 20.0,
        faces: vec![
            FaceCondition::Fixed {
                t0_c: 40.0,
                gradient: g,
            };
            6
        ],
        source_w: 0.0,
        body_sources: Vec::new(),
        hdt_c: None,
    };
    let parent =
        solve_conduction(&model, &problem, &options(1.0)).unwrap_or_else(|e| panic!("{e}"));
    let sub = solve_thermal_submodel(
        &model,
        &problem,
        &parent,
        SolidId(0),
        ([-1.0, 1.1, -1.0], [3.3, 5.0, 4.0]),
        &options(0.25),
    )
    .unwrap_or_else(|e| panic!("{e}"));
    // The x upper and y lower cuts are capped.
    assert_eq!(sub.caps.len(), 2, "{:?}", sub.caps);
    let exact = |p: &[f64; 3]| 40.0 + (0..3).map(|a| g[a] * p[a]).sum::<f64>();
    let worst = sub.model.solids[0]
        .vertices
        .iter()
        .zip(&sub.result.vertex_temperature[0])
        .map(|(p, t)| (t - exact(p)).abs())
        .fold(0.0, f64::max);
    assert!(worst <= 1e-8 * 40.0, "{worst}");
    // The flux is exact and the heat through the cut balances.
    let q = 1e3 * 20.0 * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    let flux = sub.result.vertex_heat_flux[0]
        .iter()
        .map(|v| ((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() - q).abs())
        .fold(0.0, f64::max);
    assert!(flux <= 1e-8 * q, "{flux} of {q}");
    let heat = &sub.result.heat;
    assert!(heat.imbalance().abs() <= 1e-9 * heat.throughput());
}

/// T3's tube (bore radius 2, wall 4, 100 °C inside, 20 °C outside):
/// the coarse solve's (4 cells through the wall) and a submodel's (8
/// cells per bore radius, its box reaching `half` either side of a
/// point on the bore) radial flux at the bore vertices clear of the
/// caps, over the exact k ΔT / (r ln(ro/ri)): (lowest, highest).
fn bore_flux(half: f64) -> ((f64, f64), (f64, f64)) {
    let (ri, ro, ti, to, k) = (2.0, 6.0, 100.0, 20.0, 30.0);
    let model = model_of(tube(ri, ro, 4.0, 96), 3);
    let mut faces = vec![FaceCondition::Insulated; 3];
    for (f, t) in [(0, ti), (1, to)] {
        faces[f] = FaceCondition::Fixed {
            t0_c: t,
            gradient: [0.0; 3],
        };
    }
    let problem = ConductionProblem {
        conductivity: k,
        faces,
        source_w: 0.0,
        body_sources: Vec::new(),
        hdt_c: None,
    };
    let parent =
        solve_conduction(&model, &problem, &options(1.0)).unwrap_or_else(|e| panic!("{e}"));
    let wall = 1e3 * k * (ti - to) / (ri * (ro / ri).ln());
    let h = ri / 8.0;
    let sub = solve_thermal_submodel(
        &model,
        &problem,
        &parent,
        SolidId(0),
        ([ri - half, -half, -1.0], [ri + half, half, 5.0]),
        &options(h),
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let radial = |vertices: &[[f64; 3]], flux: &[[f64; 3]], margin: f64, caps: bool| {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for (v, q) in vertices.iter().zip(flux) {
            let rad = v[0].hypot(v[1]);
            let clear = !caps
                || crate::model::cut_model::clear_of_caps(&sub.caps, sub.lo, sub.hi, *v, margin);
            if (rad - ri).abs() < 1e-6 && clear {
                let r = (q[0] * v[0] + q[1] * v[1]) / rad / wall;
                lo = lo.min(r);
                hi = hi.max(r);
            }
        }
        (lo, hi)
    };
    let coarse = radial(
        &model.solids[0].vertices,
        &parent.vertex_heat_flux[0],
        0.0,
        false,
    );
    let fine = radial(
        &sub.model.solids[0].vertices,
        &sub.result.vertex_heat_flux[0],
        crate::structural::submodel::CAP_CLEARANCE_CELLS * h,
        true,
    );
    (coarse, fine)
}

/// Measured 2026-10-11, radial flux over the exact at the bore: the
/// coarse solve 0.768..0.872; submodels 0.950..0.986 with the box 1.25
/// radii out, 0.966..0.990 at 1.75 and 0.974..0.999 at 2.5; the whole
/// tube at the submodel's cell size 0.978..0.993.
#[test]
fn a_submodel_recovers_the_bore_wall_flux() {
    let ((clo, chi), (near_lo, _)) = bore_flux(2.5);
    let (_, (lo, hi)) = bore_flux(5.0);
    assert!(clo < 0.8 && chi < 0.9, "coarse {clo}..{chi}");
    assert!(lo >= 0.97 && hi <= 1.03, "submodel {lo}..{hi}");
    // A bigger box only moves the caps' error further off.
    assert!(lo > near_lo, "{near_lo} -> {lo}");
}

/// PH-036: T3's bore through `refine_flux` (the `--refine` path, its
/// default box and resolution) at 4 cells through the wall, over four
/// grid shifts: the peak flux over k ΔT / (ri ln(ro/ri)). Measured
/// 2026-10-11: coarse 0.860–0.882, refined 0.993 at every shift (the
/// box, 4 radii out, holds the whole tube).
#[test]
fn t3_the_bore_flux_refines_within_3_percent() {
    let (ri, ro, ti, to, k) = (2.0, 6.0, 100.0, 20.0, 30.0);
    let model = model_of(tube(ri, ro, 4.0, 96), 3);
    let mut faces = vec![FaceCondition::Insulated; 3];
    for (f, t) in [(0, ti), (1, to)] {
        faces[f] = FaceCondition::Fixed {
            t0_c: t,
            gradient: [0.0; 3],
        };
    }
    let problem = ConductionProblem {
        conductivity: k,
        faces,
        source_w: 0.0,
        body_sources: Vec::new(),
        hdt_c: None,
    };
    let wall = 1e3 * k * (ti - to) / (ri * (ro / ri).ln());
    let ratios: Vec<(f64, f64)> = [
        [0.0; 3],
        [0.31, 0.17, 0.0],
        [0.53, 0.71, 0.0],
        [0.11, 0.43, 0.0],
    ]
    .iter()
    .map(|s| {
        let o = AnalysisOptions {
            grid_shift: *s,
            ..options(1.0)
        };
        let parent = solve_conduction(&model, &problem, &o).unwrap_or_else(|e| panic!("{e}"));
        let r = refine_flux(&model, &problem, &parent, (FaceId(0), ri), &o)
            .unwrap_or_else(|e| panic!("{e}"));
        (r.coarse_peak / wall, r.peak_heat_flux / wall)
    })
    .collect();
    assert!(
        ratios.iter().all(|x| (x.1 - 1.0).abs() <= 0.03),
        "{ratios:?}"
    );
}
