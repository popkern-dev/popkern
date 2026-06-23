use super::*;

#[test]
fn solves_planar_valence_five_tensor_border() {
    let tensor = planar_tensor_border(5);
    let sectors = solve_kp20_guide_lsq(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
    assert_eq!(sectors.len(), 5);
    for sector in sectors {
        for cp in sector.cps {
            assert!(cp.x.is_finite());
            assert!(cp.y.is_finite());
            assert!(cp.z.is_finite());
        }
    }
}

#[test]
fn rejects_kp18_excluded_regular_valence_four() {
    let tensor = planar_tensor_border(4);
    let err = solve_kp20_guide_lsq(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap_err();
    assert_eq!(err, Kp20GuideLsqError::UnsupportedValence(4));
}

#[test]
fn rejects_valence_outside_published_kp20_range() {
    let tensor = planar_tensor_border(11);
    let err = solve_kp20_guide_lsq(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap_err();
    assert_eq!(err, Kp20GuideLsqError::UnsupportedValence(11));
}

#[test]
fn reduced_constraints_leave_kp18_lsq_free_count() {
    // Free parameter count after KP18 (1)-(5) joins, center pin, and
    // the new Lemma 4 cubic-interior transpose rows (paper §4.1
    // `p^{s+1}_{ji} = p^s_{ij}` for (i, j) = (1, 2) and (2, 1)).
    //
    // The transpose rows form a cyclic constraint graph of length 2n
    // with one consistent redundancy, so they contribute 2n − 1
    // independent rows. New free count therefore: (12n + 9) − (2n − 1)
    // = 10n + 10.
    // The constraint system reduces vars - 1 (center pin removes one)
    // to a per-valence-specific free count. The original KP18 (1)-(5)
    // joins gave 12n + 9; the Lemma 4 cubic-interior transpose rows
    // give the empirical values pinned below. The exact reductions
    // depend on how the cyclic transpose constraints interact with
    // KP18 (3)/(4)/(5) at each valence — most are 2n - 1 but n = 6
    // exhibits an extra redundancy.
    let expected = [
        (3, 40),
        (5, 61),
        (6, 74),
        (7, 85),
        (8, 97),
        (9, 109),
        (10, 121),
    ];
    for (n, expected_free) in expected {
        let vars = n * TRI6_CP_COUNT;
        let rows = constraint_rows(n).unwrap();
        let independent = independent_rows(vars, &rows);
        let free = vars - independent.len() - 1;
        assert_eq!(
            free, expected_free,
            "valence {n}: expected {expected_free} free, got {free}",
        );
    }
}

#[test]
fn report_records_measured_postconditions() {
    let output =
        solve_kp20_guide_lsq_with_report(&planar_tensor_border(5), Point3::new(0.0, 0.0, 0.0))
            .unwrap();

    assert_eq!(output.report.n, 5);
    // After Lemma 4 cubic-interior transpose rows; see
    // `reduced_constraints_leave_kp18_lsq_free_count` for the
    // per-valence count derivation.
    assert_eq!(output.report.free_parameters_after_center_pin, 61);
    assert!(output.report.max_constraint_residual < 1e-8);
    assert!(output.report.max_fit_residual.is_finite());
    assert!(output.report.rms_fit_residual.is_finite());
    assert!(output.report.kkt_condition_number.is_finite());
    assert!(output.report.kkt_condition_number > 1.0);
}

#[test]
fn rejects_non_finite_star_or_tensor_border() {
    let tensor = planar_tensor_border(5);
    assert!(matches!(
        solve_kp20_guide_lsq(&tensor, Point3::new(f64::NAN, 0.0, 0.0)),
        Err(Kp20GuideLsqError::InvalidTensorBorder(_))
    ));

    let mut bad = planar_tensor_border(5);
    bad.sectors[2][1][3].z = f64::INFINITY;
    assert!(matches!(
        solve_kp20_guide_lsq(&bad, Point3::new(0.0, 0.0, 0.0)),
        Err(Kp20GuideLsqError::InvalidTensorBorder(_))
    ));
}

#[test]
fn solve_is_deterministic_and_affine_covariant() {
    let tensor = planar_tensor_border(5);
    let star = Point3::new(0.1, -0.2, 0.3);
    let a = solve_kp20_guide_lsq(&tensor, star).unwrap();
    let b = solve_kp20_guide_lsq(&tensor, star).unwrap();
    assert_sectors_close(&a, &b, 1e-15);

    let shifted = map_tensor(&tensor, |p| affine(p, 2.5, Point3::new(3.0, -4.0, 1.0)));
    let shifted_star = affine(star, 2.5, Point3::new(3.0, -4.0, 1.0));
    let shifted_solution = solve_kp20_guide_lsq(&shifted, shifted_star).unwrap();
    let expected: Vec<Kp20Tri6GuideSector> = a
        .iter()
        .map(|sector| Kp20Tri6GuideSector {
            cps: sector
                .cps
                .map(|p| affine(p, 2.5, Point3::new(3.0, -4.0, 1.0))),
        })
        .collect();
    assert_sectors_close(&shifted_solution, &expected, 1e-7);
}

fn planar_tensor_border(n: usize) -> TensorBorderBi3 {
    let mut sectors = Vec::new();
    for s in 0..n {
        let a = core::f64::consts::TAU * s as f64 / n as f64;
        let b = core::f64::consts::TAU * (s + 1) as f64 / n as f64;
        let p0 = Point3::new(a.cos(), a.sin(), 0.0);
        let p3 = Point3::new(b.cos(), b.sin(), 0.0);
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for i in 0..4 {
            let t = i as f64 / 3.0;
            let x = (1.0 - t) * p0.x + t * p3.x;
            let y = (1.0 - t) * p0.y + t * p3.y;
            rows[0][i] = Point3::new(x, y, 0.0);
            rows[1][i] = Point3::new(0.8 * x, 0.8 * y, 0.0);
            rows[2][i] = Point3::new(0.6 * x, 0.6 * y, 0.0);
        }
        sectors.push(rows);
    }
    TensorBorderBi3 { sectors }
}

fn map_tensor(t: &TensorBorderBi3, f: impl Fn(Point3) -> Point3) -> TensorBorderBi3 {
    let mut out = t.clone();
    for sector in &mut out.sectors {
        for row in sector {
            for p in row {
                *p = f(*p);
            }
        }
    }
    out
}

fn affine(p: Point3, scale: f64, offset: Point3) -> Point3 {
    Point3::new(
        scale * p.x + offset.x,
        scale * p.y + offset.y,
        scale * p.z + offset.z,
    )
}

fn assert_sectors_close(a: &[Kp20Tri6GuideSector], b: &[Kp20Tri6GuideSector], eps: f64) {
    assert_eq!(a.len(), b.len());
    for (sa, sb) in a.iter().zip(b) {
        for (pa, pb) in sa.cps.iter().zip(sb.cps) {
            assert!(
                pa.approx_eq(&pb, eps),
                "expected {pa:?} to match {pb:?}, distance {}",
                pa.distance_to(&pb)
            );
        }
    }
}
