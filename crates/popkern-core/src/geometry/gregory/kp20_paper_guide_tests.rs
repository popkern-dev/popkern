use super::*;
use crate::geometry::gregory::tri6_flat_index;

/// Identity gate: feeding the planar characteristic border itself must
/// reproduce the identity embedding exactly — `g ∘ χ̃ = χ̃` with
/// `g = id`, which satisfies every KP18 join row, so the constrained
/// LSQ minimum is zero and unique up to the fit's null directions.
#[test]
fn characteristic_border_input_yields_the_identity_guide() {
    for n in [3, 5, 6, 8] {
        let chi = cc_characteristic_tensor_border(n).expect("chi");
        let border = Kp20PaperGuideBorder::from_characteristic(&chi);
        let solved = solve_kp20_paper_guide(&border, Point3::new(0.0, 0.0, 0.0)).expect("solve");
        // n = 3's corner sits at local (~1.35, ~1.35), so the basis
        // extrapolation amplifies roundoff; 1e-7 remains decisive
        // against the legacy surrogate's O(1e-1) misfit.
        assert!(
            solved.max_fit_residual < 1.0e-7,
            "n={n}: fit residual {:e}",
            solved.max_fit_residual,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-8,
            "n={n}: constraint residual {:e}",
            solved.max_constraint_residual,
        );
        let alpha = core::f64::consts::TAU / n as f64;
        let mut worst = 0.0_f64;
        for (k, sector) in solved.sectors.iter().enumerate() {
            let a0 = k as f64 * alpha;
            let a1 = (k + 1) as f64 * alpha;
            let spoke0 = [a0.cos(), a0.sin()];
            let spoke1 = [a1.cos(), a1.sin()];
            for i in 0..=6_usize {
                for j in 0..=(6 - i) {
                    let cp = sector.cps[tri6_flat_index(6, i, j)];
                    // Guide chart: i (u) along spoke k, j (v) along
                    // spoke k+1.
                    let u = i as f64 / 6.0;
                    let v = j as f64 / 6.0;
                    let expected = [u * spoke0[0] + v * spoke1[0], u * spoke0[1] + v * spoke1[1]];
                    worst = worst
                        .max((cp.x - expected[0]).abs())
                        .max((cp.y - expected[1]).abs())
                        .max(cp.z.abs());
                }
            }
        }
        eprintln!(
            "identity gate n={n}: fit={fit:e} constraints={con:e} worst_cp_delta={worst:e} kkt={kkt:e}",
            fit = solved.max_fit_residual,
            con = solved.max_constraint_residual,
            kkt = solved.kkt_condition_number,
        );
        assert!(
            worst < 1.0e-6,
            "n={n}: solved guide deviates from the identity embedding by {worst:e}",
        );
    }
}

/// Diagnostic: which legacy constraint rows does the identity guide
/// violate? The identity is a valid C∞ rotationally-symmetric guide,
/// so any violated row over-constrains beyond KP18.
#[test]
fn identity_guide_constraint_audit_is_reported() {
    use crate::geometry::gregory::kp20_bi6_guide_lsq::var_index;
    let n = 5;
    let alpha = core::f64::consts::TAU / n as f64;
    let mut identity = vec![[0.0_f64; 2]; n * 28];
    for k in 0..n {
        let a0 = k as f64 * alpha;
        let a1 = (k + 1) as f64 * alpha;
        let spoke0 = [a0.cos(), a0.sin()];
        let spoke1 = [a1.cos(), a1.sin()];
        for i in 0..=6_usize {
            for j in 0..=(6 - i) {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                identity[var_index(k, tri6_flat_index(6, i, j))] =
                    [u * spoke0[0] + v * spoke1[0], u * spoke0[1] + v * spoke1[1]];
            }
        }
    }
    let rows = kp18_guide_constraint_rows(n);
    let mut violated = 0;
    let mut worst = 0.0_f64;
    let mut worst_row = 0_usize;
    for (index, row) in rows.iter().enumerate() {
        let mut acc = [0.0_f64; 2];
        for &(var, w) in &row.terms {
            acc[0] += w * identity[var][0];
            acc[1] += w * identity[var][1];
        }
        let r = acc[0].abs().max(acc[1].abs());
        if r > 1.0e-9 {
            violated += 1;
        }
        if r > worst {
            worst = r;
            worst_row = index;
        }
    }
    eprintln!(
        "identity constraint audit n={n}: violated={violated}/{total} worst={worst:e} at row {worst_row}",
        total = rows.len(),
    );
    assert_eq!(
        violated, 0,
        "the identity guide must satisfy every KP18 constraint row (worst {worst:e})",
    );
}

/// Diagnostic: evaluate the fit rows at the exact identity guide and
/// locate the worst row.
#[test]
fn identity_guide_fit_row_audit_is_reported() {
    let n = 5;
    let chi = cc_characteristic_tensor_border(n).expect("chi");
    let border = Kp20PaperGuideBorder::from_characteristic(&chi);
    let alpha = core::f64::consts::TAU / n as f64;
    // Identity guide coefficients per sector.
    let mut identity = vec![Point3::new(0.0, 0.0, 0.0); n * 28];
    for k in 0..n {
        let a0 = k as f64 * alpha;
        let a1 = (k + 1) as f64 * alpha;
        let spoke0 = [a0.cos(), a0.sin()];
        let spoke1 = [a1.cos(), a1.sin()];
        for i in 0..=6_usize {
            for j in 0..=(6 - i) {
                let u = i as f64 / 6.0;
                let v = j as f64 / 6.0;
                identity[crate::geometry::gregory::kp20_bi6_guide_lsq::var_index(
                    k,
                    tri6_flat_index(6, i, j),
                )] = Point3::new(
                    u * spoke0[0] + v * spoke1[0],
                    u * spoke0[1] + v * spoke1[1],
                    0.0,
                );
            }
        }
    }
    // Rebuild the fit rows exactly as the solver does.
    let (rows, targets) = kp20_paper_guide_fit_rows_for_test(&border);
    let mut worst = 0.0_f64;
    let mut worst_index = 0_usize;
    for (index, row) in rows.iter().enumerate() {
        let mut acc = [0.0_f64; 2];
        for &(var, w) in &row.terms {
            acc[0] += w * identity[var].x;
            acc[1] += w * identity[var].y;
        }
        let r = (acc[0] - targets[0][index])
            .abs()
            .max((acc[1] - targets[1][index]).abs());
        if r > worst {
            worst = r;
            worst_index = index;
        }
    }
    // Row layout: sector-major, then piece (2), end (2), depth (3), along (3).
    let per_sector = 2 * 2 * 9;
    let sector = worst_index / per_sector;
    let rem = worst_index % per_sector;
    let piece = rem / 18;
    let rem = rem % 18;
    let end = rem / 9;
    let rem = rem % 9;
    let depth = rem / 3;
    let along = rem % 3;
    eprintln!(
        "identity fit audit n={n}: worst={worst:e} at sector={sector} piece={piece} end={end} depth={depth} along={along}"
    );
}

/// Linear-lift gate: a plane `z = αx + βy` over the characteristic
/// border is exactly representable (affine in the sector chart), so
/// the fit must be exact and the solved guide must evaluate onto the
/// plane.
#[test]
fn linear_lift_input_is_fit_exactly() {
    let n = 5;
    let (alpha_lift, beta_lift) = (0.375, -0.2);
    let chi = cc_characteristic_tensor_border(n).expect("chi");
    let mut border = Kp20PaperGuideBorder::from_characteristic(&chi);
    for sector in &mut border.sectors {
        for piece in sector.iter_mut() {
            for row in piece.rows.iter_mut() {
                for cp in row.iter_mut() {
                    *cp = Point3::new(cp.x, cp.y, alpha_lift * cp.x + beta_lift * cp.y);
                }
            }
        }
    }
    let solved = solve_kp20_paper_guide(&border, Point3::new(0.0, 0.0, 0.0)).expect("solve");
    assert!(
        solved.max_fit_residual < 1.0e-8,
        "fit residual {:e}",
        solved.max_fit_residual,
    );
    // The solved guide must lie on the plane wherever evaluated.
    let mut worst = 0.0_f64;
    for sector in &solved.sectors {
        for &(a, b) in &[(0.2, 0.3), (0.5, 0.1), (0.7, 0.7), (0.1, 0.8)] {
            let value = sector.jet3(a, b).h[0][0];
            worst = worst.max((value.z - alpha_lift * value.x - beta_lift * value.y).abs());
        }
    }
    assert!(worst < 1.0e-7, "solved guide leaves the plane by {worst:e}");
}

/// Production-path smoke: a smooth single-cubic-per-sector border is
/// split and fit; constraints must close, the fit residual is a
/// recorded measurement of how well g∘χ̃ can represent a
/// non-characteristic border.
#[test]
fn production_border_fit_is_measured() {
    use crate::geometry::gregory::kp_polar_boundary_input::TensorBorderBi3;
    for n in [5, 6] {
        let mut sectors = Vec::with_capacity(n);
        for k in 0..n {
            let a0 = core::f64::consts::TAU * k as f64 / n as f64;
            let a1 = core::f64::consts::TAU * (k + 1) as f64 / n as f64;
            let p0 = [a0.cos(), a0.sin()];
            let p1 = [a1.cos(), a1.sin()];
            let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for (row_index, row) in rows.iter_mut().enumerate() {
                let radius = 1.0 - 0.1 * row_index as f64;
                let height = 0.15 * row_index as f64;
                for (i, cp) in row.iter_mut().enumerate() {
                    let t = i as f64 / 3.0;
                    *cp = Point3::new(
                        radius * ((1.0 - t) * p0[0] + t * p1[0]),
                        radius * ((1.0 - t) * p0[1] + t * p1[1]),
                        height,
                    );
                }
            }
            sectors.push(rows);
        }
        let border = Kp20PaperGuideBorder::from_tensor_border(&TensorBorderBi3 { sectors });
        let solved = solve_kp20_paper_guide(&border, Point3::new(0.0, 0.0, 0.45)).expect("solve");
        eprintln!(
            "production fit n={n}: fit={fit:e} constraints={con:e} kkt={kkt:e}",
            fit = solved.max_fit_residual,
            con = solved.max_constraint_residual,
            kkt = solved.kkt_condition_number,
        );
        assert!(
            solved.max_constraint_residual < 1.0e-8,
            "n={n}: constraint residual {:e}",
            solved.max_constraint_residual,
        );
        assert!(
            solved.max_fit_residual < 0.5,
            "n={n}: fit residual {:e}",
            solved.max_fit_residual,
        );
    }
}
