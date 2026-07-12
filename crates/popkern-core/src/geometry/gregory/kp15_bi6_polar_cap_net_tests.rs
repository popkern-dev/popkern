//! PK-087 stage 3 — leg-1 acceptance for the shared-net polar bridge.
//!
//! Takes a **real** polar fixture (not the synthetic oracle net),
//! runs it through [`convert_polar_to_bi6_cap_net`] and the verified
//! populator [`build_bi6_cap_from_net`], and asserts the three
//! bridge-acceptance properties:
//!
//! - **(a) no fold** — every assembled bi-9 sub-patch clears the
//!   Jacobian/self-intersection detector used by the F₅ ladder.
//! - **(b) G¹/G² at the oracle's closure level** — the cross-sector
//!   eq (2)/(3) residuals, measured in the *script convention* with
//!   the oracle's own `sector_residuals` formula, close to ≤ 1e-9.
//! - **(c) eq (2) does not blow up** — the same eq (2) residual is
//!   the historical diagnostic; the pre-existing fold path gave
//!   ~10–11 and the bad rim-as-spoke wire gave ~73.

use super::super::kp15_bi6_deep_interior::{
    Bi6DeepInteriorFree, assemble_bi6_polar_cap_from_net, build_bi6_cap_from_net,
};
use super::super::kp15_bi6_reparam::Bi6Reparams;
use super::super::rectangular_binonic::RectangularGregoryBiNonic;
use super::super::{
    GregoryPatch, OuterCornerData, PolarBoundaryInput, SpokeBoundaryData, TangentPlane,
};
use super::*;
use crate::geometry::{Point3, Vec3};

const SELF_INT_TOL: f64 = 1e-6;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn tp_xy() -> TangentPlane {
    TangentPlane::from_basis(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

/// A real polar fixture: an `n`-gon vertex blend. Spokes are
/// degree-7 curves from the centre out to unit-circle corners,
/// bowed up in +z toward the apex so the cap region is a shallow
/// dome (a non-degenerate, non-planar, geometrically honest input).
/// `apex_z` sets the central rise; `radii` set per-spoke fillet radius.
fn dome_polar_input(n: usize, apex_z: f64, radii: &[f64]) -> PolarBoundaryInput {
    let center = pt(0.0, 0.0, apex_z);
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut face_tangent_planes = Vec::with_capacity(n);
    let mut per_spoke_radii = Vec::with_capacity(n);
    for k in 0..n {
        let theta = std::f64::consts::TAU * (k as f64) / (n as f64);
        let v = pt(theta.cos(), theta.sin(), 0.0);
        // Degree-7 spoke: interpolate centre → corner, with a gentle
        // parabolic z bow so the dome is smooth (not a cone).
        let mut boundary_cps = [pt(0.0, 0.0, 0.0); 8];
        for (i, cp) in boundary_cps.iter_mut().enumerate() {
            let t = i as f64 / 7.0;
            let z = apex_z * (1.0 - t) * (1.0 + t); // apex_z at t=0, 0 at t=1
            *cp = pt(t * v.x, t * v.y, z);
        }
        let r = radii[k];
        spokes.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_c: Vec3::new(0.0, -1.0, 0.0),
            xtan_a_at_v: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_v: Vec3::new(0.0, -1.0, 0.0),
            radius_at_c: r,
            radius_at_v: r,
            edge_id: format!("edge_{k}"),
        });
        outer_corners.push(OuterCornerData {
            position: v,
            tangent_plane: tp_xy(),
        });
        face_tangent_planes.push(tp_xy());
        per_spoke_radii.push(r);
    }
    PolarBoundaryInput {
        central_vertex: center,
        central_tangent_plane: tp_xy(),
        spokes,
        outer_corners,
        face_tangent_planes,
        per_spoke_radii,
    }
}

// ---- Oracle-convention G1/G2 residual (script convention) ----
//
// Port of `n6_oracle.py::sector_residuals`. Grids are SCRIPT
// convention (extraordinary point at [6][6]). For the cross-sector
// edge between sector s and t = s+1:
//   F built from P[s]: cp_F(i,j) = P[s][6-j][6-i]
//   G built from P[t]: cp_G(i,j) = P[t][6-i][6-j]
// evaluated at v = 0. Only depth j ∈ {0,1,2} matters at v = 0.

fn bern_eval(coeffs: &[f64], u: f64) -> f64 {
    // de Casteljau on a Bernstein-coefficient vector.
    let mut w = coeffs.to_vec();
    let m = w.len();
    for r in 1..m {
        for i in 0..(m - r) {
            w[i] = (1.0 - u) * w[i] + u * w[i + 1];
        }
    }
    w[0]
}

/// Degree-`d` → degree-`d−1` Bernstein derivative coefficients:
/// `c'_i = d·(c_{i+1} − c_i)`.
fn bern_diff(coeffs: &[f64]) -> Vec<f64> {
    let d = (coeffs.len() - 1) as f64;
    (0..coeffs.len() - 1)
        .map(|i| d * (coeffs[i + 1] - coeffs[i]))
        .collect()
}

/// One coordinate channel of the eq (2)/(3) residual at the s|t edge,
/// sampled at `u`. `chan` picks x/y/z.
#[allow(clippy::too_many_arguments)]
fn residuals_channel(
    grids: &[[[Point3; 7]; 7]],
    s: usize,
    n: usize,
    c: f64,
    tau: f64,
    u: f64,
    chan: usize,
) -> (f64, f64) {
    let t = (s + 1) % n;
    let comp = |p: Point3| -> f64 {
        match chan {
            0 => p.x,
            1 => p.y,
            _ => p.z,
        }
    };
    let cp_f = |i: usize, j: usize| comp(grids[s][6 - j][6 - i]);
    let cp_g = |i: usize, j: usize| comp(grids[t][6 - i][6 - j]);

    // U-Bernstein coefficient vectors (degree 6) at v = 0.
    let f_row0: Vec<f64> = (0..7).map(|i| cp_f(i, 0)).collect();
    let fv: Vec<f64> = (0..7).map(|i| 6.0 * (cp_f(i, 1) - cp_f(i, 0))).collect();
    let fvv: Vec<f64> = (0..7)
        .map(|i| 30.0 * (cp_f(i, 2) - 2.0 * cp_f(i, 1) + cp_f(i, 0)))
        .collect();
    let gv: Vec<f64> = (0..7).map(|i| 6.0 * (cp_g(i, 1) - cp_g(i, 0))).collect();
    let gvv: Vec<f64> = (0..7)
        .map(|i| 30.0 * (cp_g(i, 2) - 2.0 * cp_g(i, 1) + cp_g(i, 0)))
        .collect();

    let val_fu = bern_eval(&bern_diff(&f_row0), u);
    let val_fuu = bern_eval(&bern_diff(&bern_diff(&f_row0)), u);
    let val_fv = bern_eval(&fv, u);
    let val_fuv = bern_eval(&bern_diff(&fv), u);
    let val_fvv = bern_eval(&fvv, u);
    let val_gv = bern_eval(&gv, u);
    let val_gvv = bern_eval(&gvv, u);

    let om = 1.0 - u;
    let g1 = val_fv + val_gv - 2.0 * c * om * om * val_fu;
    let a = 24.0 * c * (tau - 1.0) / (1.0 - u + tau * u);
    let g2 = -val_gvv + val_fvv - 4.0 * c * om * om * val_fuv
        + 4.0 * c * c * om.powi(4) * val_fuu
        + c * om.powi(3) * (a * u - 8.0 * c - a) * val_fu
        + a * om * om * val_fv;
    (g1, g2)
}

/// Worst |eq(2)| (G1) and |eq(3)| (G2) over all sector edges, all
/// coordinate channels, sampled at u ∈ {0.1..0.9}.
fn worst_g1_g2(grids: &[[[Point3; 7]; 7]], params: &Bi6Reparams) -> (f64, f64) {
    let n = params.n;
    let (mut wg1, mut wg2) = (0.0_f64, 0.0_f64);
    for s in 0..n {
        for chan in 0..3 {
            for step in 1..=9 {
                let u = step as f64 / 10.0;
                let (g1, g2) = residuals_channel(grids, s, n, params.c, params.tau, u, chan);
                wg1 = wg1.max(g1.abs());
                wg2 = wg2.max(g2.abs());
            }
        }
    }
    (wg1, wg2)
}

/// Assemble the bridge output into bi-9 sub-patches and report the
/// first folded index (if any).
fn fold_index(
    params: &Bi6Reparams,
    net: &[[[Point3; 4]; 4]],
    free: &Bi6DeepInteriorFree,
) -> Option<usize> {
    let out = assemble_bi6_polar_cap_from_net(params, net, free, 1e-9).expect("assemble cap");
    for (i, grid) in out.patches_bi9.into_iter().enumerate() {
        let gp = GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }));
        if gp.self_intersects_with_tolerance(SELF_INT_TOL) {
            return Some(i);
        }
    }
    None
}

// ---- Leg-1 acceptance ----

#[test]
fn shared_net_is_c2_consistent_pentagon() {
    // The shared-net layout must make the eq (8) ring residual `sal`
    // vanish — that IS the statement that the input net is one
    // C²-consistent B-spline surface.
    let polar = dome_polar_input(5, 0.25, &[0.15; 5]);
    let (params, net, free) = convert_polar_to_bi6_cap_net(&polar).expect("bridge");
    let grids = build_bi6_cap_from_net(&params, &net, &free).expect("build");
    assert!(
        grids.sal.x.hypot(grids.sal.y).hypot(grids.sal.z) < 1e-9,
        "sal (eq(8) closure residual) must vanish for a shared net: {:?}",
        grids.sal
    );
}

#[test]
fn leg1_pentagon_uniform_closes_and_does_not_fold() {
    let polar = dome_polar_input(5, 0.25, &[0.15; 5]);
    let (params, net, free) = convert_polar_to_bi6_cap_net(&polar).expect("bridge");
    let grids = build_bi6_cap_from_net(&params, &net, &free).expect("build");

    let (g1, g2) = worst_g1_g2(&grids.sectors, &params);
    println!(
        "pentagon5 uniform: eq(2)|g1|={g1:.3e} eq(3)|g2|={g2:.3e} sal={:?}",
        grids.sal
    );

    // (b) G1/G2 at the oracle's closure level.
    assert!(g1 <= 1e-9, "eq(2)/G1 residual {g1:.3e} exceeds 1e-9");
    assert!(g2 <= 1e-9, "eq(3)/G2 residual {g2:.3e} exceeds 1e-9");
    // (c) eq(2) did NOT blow up (historical bad wire ~73, fold path ~10-11).
    assert!(
        g1 < 1.0,
        "eq(2) residual {g1:.3e} regressed toward the bad-wire regime"
    );

    // (a) no fold.
    let fold = fold_index(&params, &net, &free);
    assert!(fold.is_none(), "cap folded at sub-patch {fold:?}");
}

#[test]
fn leg1_hexagon_variable_closes_and_does_not_fold() {
    // Mixed-radius hexagon — the PK-087 acceptance signature class.
    let radii = [0.30, 0.10, 0.25, 0.12, 0.28, 0.11];
    let polar = dome_polar_input(6, 0.3, &radii);
    let (params, net, free) = convert_polar_to_bi6_cap_net(&polar).expect("bridge");
    let grids = build_bi6_cap_from_net(&params, &net, &free).expect("build");

    let (g1, g2) = worst_g1_g2(&grids.sectors, &params);
    println!(
        "hexagon6 variable: eq(2)|g1|={g1:.3e} eq(3)|g2|={g2:.3e} sal={:?}",
        grids.sal
    );

    assert!(g1 <= 1e-9, "eq(2)/G1 residual {g1:.3e} exceeds 1e-9");
    assert!(g2 <= 1e-9, "eq(3)/G2 residual {g2:.3e} exceeds 1e-9");
    assert!(
        g1 < 1.0,
        "eq(2) residual {g1:.3e} regressed toward the bad-wire regime"
    );

    let fold = fold_index(&params, &net, &free);
    assert!(fold.is_none(), "cap folded at sub-patch {fold:?}");
}

#[test]
fn residual_metric_detects_non_closure() {
    // Negative control: proves `worst_g1_g2` is measuring real
    // continuity (not trivially returning ~0). Perturbing a
    // near-spoke CP of one sector breaks the cross-sector eq (2)/(3)
    // constraints, so the residual must jump well above the closure
    // floor.
    let polar = dome_polar_input(5, 0.25, &[0.15; 5]);
    let (params, net, free) = convert_polar_to_bi6_cap_net(&polar).expect("bridge");
    let mut grids = build_bi6_cap_from_net(&params, &net, &free).expect("build");
    let (g1_clean, _) = worst_g1_g2(&grids.sectors, &params);
    assert!(g1_clean <= 1e-9, "clean baseline {g1_clean:.3e}");

    // Corrupt a first-cross-derivative CP (script [6][5], on the
    // spoke, depth 1) of sector 0 — directly in the eq (2) support.
    grids.sectors[0][6][5].x += 0.5;
    let (g1_dirty, _) = worst_g1_g2(&grids.sectors, &params);
    assert!(
        g1_dirty > 1e-3,
        "perturbing a spoke CP must break eq(2): got {g1_dirty:.3e} (baseline {g1_clean:.3e})"
    );
}

// ---- PK-087 stage 4: measured boundary-vs-continuity trade-off ----
//
// The acceptance gate is `cap_outer_g0 ≤ 5 mm` (boundary fidelity of
// the cap's outer edge to the surrounding surface). Two source
// mechanisms exist for the cap's tensor-border rows; this test pins
// the measured trade-off between them so a future C²-consistent
// net-fit (the commercial-grade replacement) has a numeric baseline.
//
//  - **Shared B-spline cap net** ([`build_bi6_cap_from_net`]): the net
//    sharing makes the input one C²-consistent surface, so eq (2)/(3)
//    close at machine zero — BUT the B-spline → Bézier averaging pulls
//    the outer boundary *inward*, so it cannot interpolate the real
//    surrounding-surface boundary. Boundary error ≫ 5 mm.
//  - **Real tensor border direct** ([`build_bi6_cap_from_tensor_border`],
//    the paper `TBdata` path): the outer boundary interpolates the
//    real boundary's endpoints exactly (error 0), and `sal` still
//    vanishes — BUT the crude per-sector fillet extractor border is
//    not C²-consistent across the shared spokes, so eq (2)/(3) do NOT
//    close (G1/G2 ≫ 5e-2).
//
// Neither pure mechanism clears both gates; the commercial-grade fix
// is a C²-consistent net *fitted* to the real boundary (recorded in
// the PK-087 ADR).

fn dist(a: Point3, b: Point3) -> f64 {
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
}

#[test]
fn net_path_closes_continuity_but_shrinks_boundary() {
    use super::super::kp_polar_boundary_input::TensorBorder;
    use super::super::kp15_bi6_tensor_border_extract::extract_bi3_tensor_border_from_polar_input;
    let radii = [0.30, 0.10, 0.25, 0.12, 0.28, 0.11];
    let polar = dome_polar_input(6, 0.3, &radii);
    let TensorBorder::Bi3(tb) = extract_bi3_tensor_border_from_polar_input(&polar).expect("tb")
    else {
        panic!("bi3");
    };
    let (params, net, free) = convert_polar_to_bi6_cap_net(&polar).expect("bridge");
    let grids = build_bi6_cap_from_net(&params, &net, &free).expect("build");
    let (g1, g2) = worst_g1_g2(&grids.sectors, &params);
    // Continuity closes at machine zero.
    assert!(g1 <= 1e-9 && g2 <= 1e-9, "net path g1={g1:.3e} g2={g2:.3e}");
    // But the outer corner (script [0][0]) is pulled well inside the
    // real corner — boundary fidelity is far above the 5 mm gate.
    let mut worst = 0.0_f64;
    for s in 0..6 {
        worst = worst.max(dist(grids.sectors[s][0][0], tb.sectors[s][0][0]));
    }
    assert!(
        worst > 0.1,
        "net path is expected to shrink the boundary (measured on this fixture ~0.52); got {worst:.3e}"
    );
}

#[test]
fn tensor_border_direct_matches_boundary_but_breaks_continuity() {
    use super::super::kp_polar_boundary_input::TensorBorder;
    use super::super::kp15_bi6_deep_interior::{
        SectorTensorBorderBi3, build_bi6_cap_from_tensor_border,
    };
    use super::super::kp15_bi6_tensor_border_extract::extract_bi3_tensor_border_from_polar_input;
    let radii = [0.30, 0.10, 0.25, 0.12, 0.28, 0.11];
    let polar = dome_polar_input(6, 0.3, &radii);
    let n = 6usize;
    let TensorBorder::Bi3(tb) = extract_bi3_tensor_border_from_polar_input(&polar).expect("tb")
    else {
        panic!("bi3");
    };
    let (params, _net, free) = convert_polar_to_bi6_cap_net(&polar).expect("bridge");
    // Shared-spoke-consistent mapping: tak(s) shares with tad(s-1).
    let tad_src: Vec<SectorTensorBorderBi3> = (0..n).map(|s| tb.sectors[s]).collect();
    let tak_src: Vec<SectorTensorBorderBi3> = (0..n).map(|s| tb.sectors[(s + n - 1) % n]).collect();
    let grids =
        build_bi6_cap_from_tensor_border(&params, &tad_src, &tak_src, &free).expect("build");
    // sal still vanishes (eq(8) ring closes) ...
    let sal = grids.sal.x.hypot(grids.sal.y).hypot(grids.sal.z);
    assert!(sal < 1e-9, "sal {sal:.3e}");
    // ... and the outer boundary interpolates the real border endpoints.
    let mut worst_border = 0.0_f64;
    for s in 0..n {
        worst_border = worst_border.max(dist(grids.sectors[s][0][0], tb.sectors[s][0][0]));
        worst_border = worst_border.max(dist(grids.sectors[s][6][0], tb.sectors[s][0][3]));
    }
    assert!(
        worst_border < 1e-12,
        "TBdata-direct boundary must interpolate the input endpoints; got {worst_border:.3e}"
    );
    // But eq(2)/(3) do NOT close — the extractor border is not
    // C²-consistent across the shared spokes.
    let (g1, _g2) = worst_g1_g2(&grids.sectors, &params);
    assert!(
        g1 > 5e-2,
        "TBdata-direct is expected to break cross-sector continuity (measured ~1.97); got g1={g1:.3e}"
    );
}

#[test]
fn declines_low_valence() {
    let polar = dome_polar_input(4, 0.2, &[0.15; 4]);
    // n = 4 has fewer than 5 spokes → shared-net bridge declines.
    // (dome_polar_input(4) is a valid PolarBoundaryInput; the bridge
    // gate is what rejects it.)
    assert!(convert_polar_to_bi6_cap_net(&polar).is_err());
}
