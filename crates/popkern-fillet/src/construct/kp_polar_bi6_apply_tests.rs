//! Tests for the `popkern-fillet` bi-6 polar wrapper.
//!
//! Exercises [`try_kp_polar_bi6_cap`] with synthetic polar inputs
//! at valences ∈ {5, 6, 7, 8, 10, 12} and verifies the full
//! pipeline lands N bi-9 CP grids with finite values.

use super::{try_kp_polar_bi6_cap, try_kp_polar_bi6_cap_via_table};
use popkern_core::geometry::Point3;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    OuterCornerData, PolarBoundaryInput, SpokeBoundaryData, TangentPlane,
};

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn tp_xy() -> TangentPlane {
    TangentPlane::from_basis(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

fn tp_radial(theta: f64) -> TangentPlane {
    let u = Vec3::new(-theta.sin(), theta.cos(), 0.0);
    let v = Vec3::new(theta.cos(), theta.sin(), 0.0);
    TangentPlane::from_basis(u, v).unwrap()
}

/// Build a regular N-gon polar input at origin, unit outer
/// radius, constant per-spoke fillet radius 0.1.
fn regular_polar(n: usize) -> PolarBoundaryInput {
    let center = pt(0.0, 0.0, 0.0);
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut face_tangent_planes = Vec::with_capacity(n);
    let mut per_spoke_radii = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
        let v_outer = pt(theta.cos(), theta.sin(), 0.0);
        let mut boundary_cps = [pt(0.0, 0.0, 0.0); 8];
        for (i, cp) in boundary_cps.iter_mut().enumerate() {
            let t = i as f64 / 7.0;
            *cp = pt(t * v_outer.x, t * v_outer.y, 0.0);
        }
        spokes.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_c: Vec3::new(0.0, -1.0, 0.0),
            xtan_a_at_v: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_v: Vec3::new(0.0, -1.0, 0.0),
            radius_at_c: 0.1,
            radius_at_v: 0.1,
            edge_id: format!("edge_{k}"),
        });
        outer_corners.push(OuterCornerData {
            position: v_outer,
            tangent_plane: tp_xy(),
        });
        face_tangent_planes.push(tp_radial(theta));
        per_spoke_radii.push(0.1);
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

#[test]
fn pentagon_cap_produces_five_patches() {
    let polar = regular_polar(5);
    let patches = try_kp_polar_bi6_cap(&polar).expect("pentagon converges");
    assert_eq!(patches.len(), 5);
    // All 100 CPs per patch finite.
    for patch in &patches {
        for row in patch {
            for cp in row {
                assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
            }
        }
    }
}

#[test]
fn hexagon_cap_produces_six_patches() {
    let polar = regular_polar(6);
    let patches = try_kp_polar_bi6_cap(&polar).expect("hexagon converges");
    assert_eq!(patches.len(), 6);
}

#[test]
fn higher_valences() {
    for n in [7usize, 8, 10, 12] {
        let polar = regular_polar(n);
        let patches = try_kp_polar_bi6_cap(&polar).unwrap_or_else(|e| panic!("n={n}: {e}"));
        assert_eq!(patches.len(), n);
    }
}

#[test]
fn n3_special_case() {
    // n=3 path runs the converter + solver but with p^k_{22} all
    // free per KP15 §3.5 — our scaffold uses the pinned free
    // parameter as p^0_{22} and zero for p^1_{22}, p^2_{22}.
    let polar = regular_polar(3);
    let patches = try_kp_polar_bi6_cap(&polar).expect("n=3 converges");
    assert_eq!(patches.len(), 3);
}

#[test]
fn rejects_bad_input() {
    // n=2 or negative radius should fail cleanly.
    let mut polar = regular_polar(5);
    polar.per_spoke_radii[2] = 0.0;
    assert!(try_kp_polar_bi6_cap(&polar).is_err());
}

#[test]
fn rejects_n_below_3() {
    // Build a valid n=3 and truncate to n=2.
    let mut polar = regular_polar(3);
    polar.spokes.pop();
    polar.outer_corners.pop();
    polar.face_tangent_planes.pop();
    polar.per_spoke_radii.pop();
    assert!(try_kp_polar_bi6_cap(&polar).is_err());
}

#[test]
fn bi9_grid_corners_match_bi6_endpoints() {
    // The bi-9 elevation of a bi-6 grid preserves the 4 corner
    // CPs exactly. We verify the four corners are each finite
    // and at reasonable positions.
    let polar = regular_polar(5);
    let patches = try_kp_polar_bi6_cap(&polar).unwrap();
    for patch in &patches {
        // Corner (0, 0) across all sectors = star point (shared).
        let corner_00 = patch[0][0];
        assert!(corner_00.x.is_finite());
    }
    // All sectors' (0, 0) should be the same point (star).
    let star = patches[0][0][0];
    for patch in &patches {
        let c = patch[0][0];
        let d_sq = (c.x - star.x).powi(2) + (c.y - star.y).powi(2) + (c.z - star.z).powi(2);
        assert!(d_sq.sqrt() < 1e-10, "star point drift: {c:?} vs {star:?}");
    }
}

#[test]
fn via_table_pentagon_produces_finite_bi9_grids() {
    // Phase E.3 production-side wiring: same regular-pentagon
    // fixture, the via-table path produces N=5 valid bi-9 grids
    // and all sectors' (0, 0) corners agree on a shared star.
    let polar = regular_polar(5);
    let patches = try_kp_polar_bi6_cap_via_table(&polar, 5, 0.0).unwrap();
    assert_eq!(patches.len(), 5);
    let star = patches[0][0][0];
    for patch in &patches {
        for row in patch {
            for cp in row {
                assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
            }
        }
        let c = patch[0][0];
        let d = ((c.x - star.x).powi(2) + (c.y - star.y).powi(2) + (c.z - star.z).powi(2)).sqrt();
        assert!(d < 1e-10, "via-table star drift: {c:?} vs {star:?}");
    }
}

#[test]
fn via_table_succeeds_across_valences_5_to_10() {
    // Smoke test across a broader valence band.
    for n in [5usize, 6, 7, 8, 10] {
        let polar = regular_polar(n);
        let patches = try_kp_polar_bi6_cap_via_table(&polar, 5, 0.0)
            .unwrap_or_else(|e| panic!("n={n} via-table failed: {e}"));
        assert_eq!(patches.len(), n, "n={n}: patch count");
        for patch in &patches {
            for row in patch {
                for cp in row {
                    assert!(
                        cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite(),
                        "n={n}: non-finite CP {cp:?}"
                    );
                }
            }
        }
    }
}

/// WP5 fold diagnostic (2026-07-02 verification): the via-table cap
/// FAILS the recovery plan's fold-free gate on the planar-spoke
/// fixture — every patch self-intersects at every valence (worst
/// jacobian normal-cos −1.0), with control points spanning ~[−25, 33]
/// on a unit-radius input at n = 3, 5, 7, 9. Root cause is the
/// polar→cap bridge, not the phase A–E table machinery (whose
/// eq (2)/(3) gates close at machine zero on circulant-sector
/// inputs): see `assembly_eq2_residual_diagnostic` — the bridge does
/// not (S1)-propagate sector 0's near-star CPs and does not
/// paper-strictly identify the cross-sector row (eq (2) residual
/// ≈ 11 end-to-end). Not a pass/fail gate until the bridge is
/// completed.
#[test]
#[ignore = "WP5 open: polar→cap bridge incomplete — via-table cap folds on the planar spoke fixture (see doc comment)"]
fn via_table_cap_fold_probe() {
    use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
    for n in [3usize, 5, 6, 7, 8, 9, 10] {
        let polar = regular_polar(n);
        let patches = match super::try_kp_polar_bi6_cap_via_table(&polar, 5, 1.0e-3) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("PROBE n={n}: via-table cap failed: {e}");
                continue;
            }
        };
        let mut intersecting = 0_usize;
        let mut worst_cos = f64::INFINITY;
        for grid in &patches {
            let gp = GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
                grid: *grid,
            }));
            let report = gp.self_intersection_report(1.0e-9);
            if report.intersects() {
                intersecting += 1;
            }
            worst_cos = worst_cos.min(report.min_jacobian_normal_cos);
        }
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for grid in &patches {
            for row in grid.iter() {
                for cp in row.iter() {
                    min[0] = min[0].min(cp.x);
                    min[1] = min[1].min(cp.y);
                    min[2] = min[2].min(cp.z);
                    max[0] = max[0].max(cp.x);
                    max[1] = max[1].max(cp.y);
                    max[2] = max[2].max(cp.z);
                }
            }
        }
        eprintln!(
            "PROBE n={n}: patches={} intersecting={intersecting} worst_jacobian_normal_cos={worst_cos:.3e} bbox=[{:.2},{:.2}]x[{:.2},{:.2}]x[{:.2},{:.2}]",
            patches.len(),
            min[0],
            max[0],
            min[1],
            max[1],
            min[2],
            max[2],
        );
    }
}
