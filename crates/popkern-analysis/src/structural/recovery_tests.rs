//! Stress-recovery gates (PH-026): C5 stress and ZZ1, C6 (Roark plate
//! with a hole), C7 (Kirsch) and ZZ1, TW1 and TW2 (thin walls).

use super::assemble::FaceMechanics;
use super::thin_wall_tests::rotate_vector;
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{
    extrude_faces, model_of, quarter_plate_with_hole, quarter_tube, rotated,
};
use crate::recover::tensor::Voigt;
use crate::recover::{Fallback, Recovery};
use crate::report::Warning;

const ROLLER: FaceMechanics = FaceMechanics::Support(PieceSupport::Roller);

fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-10,
        cg_max_iters: 400,
        ..AnalysisOptions::default()
    }
}

/// Polar stresses (rr, θθ, rθ, zz) at (x, y) to Cartesian Voigt.
fn polar_to_cartesian(x: f64, y: f64, rr: f64, tt: f64, rt: f64, zz: f64) -> Voigt {
    let r = (x * x + y * y).sqrt();
    let (c, s) = (x / r, y / r);
    [
        rr * c * c + tt * s * s - 2.0 * rt * s * c,
        rr * s * s + tt * c * c + 2.0 * rt * s * c,
        zz,
        0.0,
        0.0,
        (rr - tt) * s * c + rt * (c * c - s * s),
    ]
}

/// Mean relative hoop-stress error over every bore vertex of the
/// plane-strain quarter tube, and the ZZ effectivity, at cell size `h`.
fn lame_recovery(h: f64) -> (f64, f64) {
    let (ri, ro, p, e, nu) = (2.0, 6.0, 1e-3, 1.0, 0.3);
    let model = model_of(quarter_tube(ri, ro, 2.0, 48), 5);
    let faces = vec![
        FaceMechanics::Pressure(p),
        FaceMechanics::Free,
        ROLLER,
        ROLLER,
        ROLLER,
    ];
    let problem = ElasticProblem {
        youngs: e,
        poisson: nu,
        faces,
    };
    let a = p * ri * ri / (ro * ro - ri * ri);
    let exact = move |x: [f64; 3]| -> Voigt {
        let r2 = x[0] * x[0] + x[1] * x[1];
        let (rr, tt) = (a * (1.0 - ro * ro / r2), a * (1.0 + ro * ro / r2));
        polar_to_cartesian(x[0], x[1], rr, tt, 0.0, nu * (rr + tt))
    };
    let r = solve_elastic_with(
        &model,
        &problem,
        &options(h),
        Fallback::InCell,
        Some(&exact),
        crate::moments::Splitting::Exact,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let hoop_exact = a * (1.0 + ro * ro / (ri * ri));
    let errors: Vec<f64> = model.solids[0]
        .vertices
        .iter()
        .zip(&r.stresses[0].vertices)
        .filter(|(v, _)| ((v[0] * v[0] + v[1] * v[1]).sqrt() - ri).abs() < 1e-9)
        .map(|(v, s)| {
            let rad = ri;
            let (c, sn) = (v[0] / rad, v[1] / rad);
            let hoop = s.stress[0] * sn * sn + s.stress[1] * c * c - 2.0 * s.stress[5] * sn * c;
            (hoop - hoop_exact).abs() / hoop_exact
        })
        .collect();
    assert!(!errors.is_empty(), "no bore vertices sampled");
    let mean = errors.iter().sum::<f64>() / errors.len() as f64;
    let theta = r.stresses[0].global_error.unwrap_or(f64::NAN)
        / r.stresses[0].true_error.unwrap_or(f64::NAN);
    (mean, theta)
}

/// C5 stress part and ZZ1. Measured: bore hoop stress mean error 34%
/// at 4 cells through the wall and 5.4% at 8 (worst 12%): linear SPR
/// extrapolating to the bore under-reads a 1/r² peak, so the map's "5%
/// at 4 cells" is carried by PH-035. Gated here: convergence (≥ 3×
/// from 4 to 8 cells, ≤ 10% at 8) and the ZZ effectivity (measured
/// 1.33 and 1.06).
#[test]
fn c5_hoop_stress_converges_and_zz1_is_honest() {
    let (coarse, theta_coarse) = lame_recovery(1.0);
    let (fine, theta_fine) = lame_recovery(0.5);
    assert!(fine <= 0.10, "8 cells: {fine}");
    assert!(coarse / fine >= 3.0, "{coarse} -> {fine}");
    for theta in [theta_coarse, theta_fine] {
        assert!((0.6..=1.6).contains(&theta), "ZZ1 C5 effectivity {theta}");
    }
}

/// The quarter plate with a hole, 2 mm half-thickness (a real plate:
/// at h = 0.5 it is 4 cells through): (model, outer segments).
#[allow(clippy::type_complexity)]
pub(super) fn plate() -> (crate::model::AnalysisModel, Vec<(u32, [f64; 2], [f64; 2])>) {
    let (mesh, nfaces, segments) = quarter_plate_with_hole(18.0, 6.0, 2.0, 2.0, 48, 18);
    (model_of(mesh, nfaces), segments)
}

/// Faces 0 (bottom, mid-plane symmetry), 2 (y = 0) and 3 (x = 0) on
/// rollers; the hole and the top free.
pub(super) fn plate_faces(nfaces: usize) -> Vec<FaceMechanics> {
    let mut faces = vec![FaceMechanics::Free; nfaces];
    faces[0] = ROLLER;
    faces[2] = ROLLER;
    faces[3] = ROLLER;
    faces
}

/// σ_xx at the hole edge on the x = 0 symmetry plane (the peak point),
/// at the mid-plane z = 0 where the 3D concentration is at least the
/// plane-stress value.
fn mid_plane_peak(model: &crate::model::AnalysisModel, r: &StructuralResult) -> f64 {
    model.solids[0]
        .vertices
        .iter()
        .zip(&r.stresses[0].vertices)
        .find(|(v, _)| v[0].abs() < 1e-9 && (v[1] - 2.0).abs() < 1e-9 && v[2].abs() < 1e-9)
        .map_or(f64::NAN, |(_, s)| s.stress[0])
}

/// C6 at r = 4 cells. Measured (PH-026): mid-plane 0.84 / 0.98 / 0.97
/// of the plane-stress Kt σ_nom at r = 3 / 4 / 6 cells; the free surface
/// reads lower (0.76 / 0.88 / 0.85), as a thick plate's does. The map's
/// r = 3 cells is not met by linear SPR; filed as its own task.
#[test]
fn c6_plate_with_a_hole_matches_roark() {
    let (model, segments) = plate();
    let sigma = 1e-3;
    let mut faces = plate_faces(model.faces.len());
    for (f, _, n) in &segments {
        if n[0] > 0.5 {
            faces[*f as usize] = FaceMechanics::Traction([sigma, 0.0, 0.0]);
        }
    }
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let r = solve_elastic(&model, &problem, &options(0.5)).unwrap_or_else(|e| panic!("{e}"));
    // Roark 17.1 case 7a, 2r/D = 1/3: Kt = 2.3067 on the net section.
    let q: f64 = 1.0 / 3.0;
    let kt = 3.0 - 3.13 * q + 3.66 * q * q - 1.53 * q * q * q;
    let peak_exact = kt * sigma * 12.0 / 8.0;
    let peak = mid_plane_peak(&model, &r);
    assert!(
        (peak - peak_exact).abs() <= 0.10 * peak_exact,
        "peak {peak} vs {peak_exact}"
    );
}

/// Kirsch's field around a hole of radius `a` in a plate under remote
/// tension `t` along x.
pub(super) fn kirsch(t: f64, a: f64) -> impl Fn(f64, f64) -> Voigt + Copy {
    move |x: f64, y: f64| -> Voigt {
        let r2 = x * x + y * y;
        let (a2, a4) = (a * a / r2, a * a * a * a / (r2 * r2));
        let th = y.atan2(x);
        let (c2, s2) = ((2.0 * th).cos(), (2.0 * th).sin());
        let rr = 0.5 * t * (1.0 - a2) + 0.5 * t * (1.0 - 4.0 * a2 + 3.0 * a4) * c2;
        let tt = 0.5 * t * (1.0 + a2) - 0.5 * t * (1.0 + 3.0 * a4) * c2;
        let rt = -0.5 * t * (1.0 + 2.0 * a2 - 3.0 * a4) * s2;
        polar_to_cartesian(x, y, rr, tt, rt, 0.0)
    }
}

/// C7: the quarter plate loaded on its outer faces by Kirsch's field.
pub(super) fn kirsch_problem(
    model: &crate::model::AnalysisModel,
    segments: &[(u32, [f64; 2], [f64; 2])],
    t: f64,
) -> ElasticProblem {
    let kirsch = kirsch(t, 2.0);
    let mut faces = plate_faces(model.faces.len());
    for (f, mid, n) in segments {
        let s = kirsch(mid[0], mid[1]);
        let traction = [s[0] * n[0] + s[5] * n[1], s[5] * n[0] + s[1] * n[1], 0.0];
        faces[*f as usize] = FaceMechanics::Traction(traction);
    }
    ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    }
}

#[test]
fn c7_kirsch_peak_and_zz1() {
    let (model, segments) = plate();
    let t = 1e-3;
    let kirsch = kirsch(t, 2.0);
    let problem = kirsch_problem(&model, &segments, t);
    let exact = move |x: [f64; 3]| kirsch(x[0], x[1]);
    let r = solve_elastic_with(
        &model,
        &problem,
        &options(0.5),
        Fallback::InCell,
        Some(&exact),
        crate::moments::Splitting::Exact,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let peak = mid_plane_peak(&model, &r);
    assert!(
        (peak - 3.0 * t).abs() <= 0.10 * 3.0 * t,
        "Kirsch peak {peak} vs {}",
        3.0 * t
    );
    let (zz, truth) = (r.stresses[0].global_error, r.stresses[0].true_error);
    let theta = zz.unwrap_or(f64::NAN) / truth.unwrap_or(f64::NAN);
    assert!(
        (0.6..=1.6).contains(&theta),
        "ZZ1 C7 effectivity {theta} ({zz:?} / {truth:?})"
    );
}

/// A wall: x in [0, 10], y in [0, 4], z in [0, t], rotated 17°, clamped
/// at x = 0 (face 4), end shear on x = 10 (face 3). Returns the model
/// and the end force.
fn wall(t: f64) -> (crate::model::AnalysisModel, ElasticProblem, f64) {
    let n = 20;
    let mut outline: Vec<[f64; 2]> = (0..=n).map(|i| [10.0 * i as f64 / n as f64, 0.0]).collect();
    outline.extend((0..=n).rev().map(|i| [10.0 * i as f64 / n as f64, 4.0]));
    let mut edge_faces = vec![2u32; outline.len()];
    edge_faces[n] = 3;
    edge_faces[2 * n + 1] = 4;
    // Strip triangulation between the two long edges (a fan from a
    // corner would put zero-area triangles along its own edge).
    let top = |i: usize| (2 * n + 1 - i) as u32;
    let tris: Vec<[u32; 3]> = (0..n)
        .flat_map(|i| {
            [
                [i as u32, (i + 1) as u32, top(i + 1)],
                [i as u32, top(i + 1), top(i)],
            ]
        })
        .collect();
    let mesh = rotated(
        extrude_faces(&outline, &tris, 0.0, t, &edge_faces),
        17.0,
        [0.0, 2.0, 0.5 * t],
    );
    let model = model_of(mesh, 5);
    let traction = 1e-4;
    let mut faces = vec![FaceMechanics::Free; 5];
    faces[4] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[3] = FaceMechanics::Traction(rotate_vector([0.0, 0.0, -traction], 17.0));
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    (model, problem, traction * 4.0 * t)
}

/// Bending stress σ_xx (wall frame) at mid-span top-surface vertices,
/// with how each was recovered.
fn mid_span_surface(
    model: &crate::model::AnalysisModel,
    r: &StructuralResult,
    t: f64,
) -> Vec<(f64, Recovery)> {
    let centre = [0.0, 2.0, 0.5 * t];
    let ex = rotate_vector([1.0, 0.0, 0.0], 17.0);
    model.solids[0]
        .vertices
        .iter()
        .zip(&r.stresses[0].vertices)
        .filter_map(|(v, s)| {
            let local = rotate_vector(
                [v[0] - centre[0], v[1] - centre[1], v[2] - centre[2]],
                -17.0,
            );
            let at_mid = (local[0] - 5.0).abs() < 1e-6 && (local[2] - 0.5 * t).abs() < 1e-6;
            at_mid.then(|| {
                let st = s.stress;
                let m = [
                    [st[0], st[5], st[4]],
                    [st[5], st[1], st[3]],
                    [st[4], st[3], st[2]],
                ];
                let mx: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| m[i][j] * ex[j]).sum());
                ((0..3).map(|i| ex[i] * mx[i]).sum(), s.recovery)
            })
        })
        .collect()
}

/// TW1. Measured (PH-026) at mid-span of a 17°-rotated wall 1.5 cells
/// thick: production recovery 1.013 / 1.027 of beam theory (patches
/// with only full cells are accurate here), the in-cell EAS field at
/// every vertex 1.003 / 0.931, constant (centroid) values 0.827 / 0.559.
#[test]
fn tw1_a_thin_rotated_wall_in_bending_keeps_its_surface_stress() {
    let t = 1.5;
    let (model, problem, force) = wall(t);
    // σ = M c / I at mid-span, tension on top for a downward end load.
    let exact = force * 5.0 * (0.5 * t) / (4.0 * t * t * t / 12.0);
    let worst = |fallback: Fallback| -> f64 {
        let r = solve_elastic_with(
            &model,
            &problem,
            &options(1.0),
            fallback,
            None,
            crate::moments::Splitting::Exact,
        )
        .unwrap_or_else(|e| panic!("{e}"));
        let samples = mid_span_surface(&model, &r, t);
        assert!(!samples.is_empty());
        samples
            .iter()
            .map(|(s, _)| (s - exact).abs() / exact)
            .fold(0.0, f64::max)
    };
    assert!(worst(Fallback::InCell) <= 0.10, "production recovery");
    assert!(worst(Fallback::ForceInCell) <= 0.10, "the in-cell fallback");
    // Negative control: constant (centroid) values sit near the neutral
    // axis on a wall in bending.
    assert!(
        worst(Fallback::ForceConstant) > 0.10,
        "constant values passed"
    );
}

#[test]
fn tw2_walls_under_one_and_a_half_cells_are_under_resolved() {
    let thin = |t: f64| -> (bool, bool) {
        let (model, problem, _) = wall(t);
        let r = solve_elastic(&model, &problem, &options(1.0)).unwrap_or_else(|e| panic!("{e}"));
        let warned = r
            .report
            .warnings
            .iter()
            .any(|w| matches!(w, Warning::UnderResolved { .. }));
        let no_bars = r.stresses[0].vertices.iter().all(|v| v.error.is_none());
        (warned, no_bars)
    };
    assert_eq!(thin(1.4), (true, true));
    // A 0.3-cell wall warns too: the rules that keep wedges beside acute
    // edges and mismatched faces quiet (PH-055) do not hide real walls.
    assert_eq!(thin(0.3), (true, true));
    let (warned, _) = thin(3.0);
    assert!(!warned, "a 3-cell wall is resolved");
}

/// X8/E2, stress part, at a regular point (mid-span of the rotated
/// 3-cell wall; a singular point such as a clamped corner has no
/// converged peak). Measured over ten grid shifts (PH-069): 0.983..1.015
/// of beam theory, a 3.2% spread, with normal extrapolation at flat and
/// convex surfaces (interpolated nodal σ* alone: 0.953..1.035, 8.2%).
/// Gated: every value within 5% of beam theory and spread under 5%.
#[test]
fn x8_surface_stress_is_steady_over_grid_shifts() {
    let t = 3.0;
    let (model, problem, force) = wall(t);
    let exact = force * 5.0 * (0.5 * t) / (4.0 * t * t * t / 12.0);
    let mut values = Vec::new();
    for k in 0..10 {
        let f = k as f64 / 10.0;
        let o = AnalysisOptions {
            grid_shift: [0.37 * f, 0.61 * f, 0.83 * f],
            ..options(1.0)
        };
        let r = solve_elastic(&model, &problem, &o).unwrap_or_else(|e| panic!("{e}"));
        values.extend(mid_span_surface(&model, &r, t).iter().map(|(s, _)| *s));
    }
    assert_eq!(values.len(), 20, "two mid-span surface vertices per shift");
    let (lo, hi) = values
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &v| {
            (l.min(v), h.max(v))
        });
    assert!(
        values.iter().all(|v| (v - exact).abs() <= 0.05 * exact),
        "{lo}..{hi} vs {exact}"
    );
    assert!(hi - lo <= 0.05 * exact, "spread {lo}..{hi}");
}
