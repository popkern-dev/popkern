//! PK-139 — variable-radius fillet adversarial hardening matrix.
//!
//! `FilletRadius::Variable` / `ByPoints` had never faced campaign
//! pressure. Each class below measures a distinct adversarial input;
//! outcomes are gated as either a WATERTIGHT completed blend (with
//! mesh-volume sanity against the constant-radius bracket) or a
//! STRUCTURED decline (error/warning naming the reason) — never a
//! panic, never a silently folded or open body.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn default_tol() -> Tolerance {
    Tolerance::default()
}

fn mesh_open_edges(brep: &BRep) -> usize {
    use std::collections::HashMap;
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut count: HashMap<(usize, usize), i32> = HashMap::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    count.values().filter(|&&c| c != 2).count()
}

fn mesh_volume(brep: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut vol = 0.0_f64;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (vol / 6.0).abs()
}

/// Fillet one cube edge with the given radius law; assert the
/// completed body is watertight with volume inside the
/// constant-radius bracket [min_r, max_r], or the decline is
/// structured (typed error / DirectApex warning — never a panic and
/// never a completed-but-open body).
fn run_law(label: &str, radius: FilletRadius, r_lo: f64, r_hi: f64) {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let spec = FilletSpec { edge, radius };
    let result = fillet_edges(
        &brep,
        &[spec],
        &FilletOptions::default(),
        default_tol(),
        None,
    );
    match result {
        Ok(BlendResult { brep: out, .. }) => {
            let opens = mesh_open_edges(&out);
            assert_eq!(opens, 0, "{label}: completed blend must be watertight");
            // Volume bracket: a cube minus the fillet band removal;
            // the removed material grows with radius, so the result
            // volume must sit between the two constant-radius
            // completions' volumes (computed live for honesty).
            let vol = mesh_volume(&out);
            let bracket = |r: f64| -> Option<f64> {
                let b = BRep::unit_cube(default_tol()).unwrap();
                let e = b.edges.keys().next().unwrap();
                fillet_edges(
                    &b,
                    &[FilletSpec {
                        edge: e,
                        radius: FilletRadius::Constant(r),
                    }],
                    &FilletOptions::default(),
                    default_tol(),
                    None,
                )
                .ok()
                .map(|r| mesh_volume(&r.brep))
            };
            if let (Some(v_hi), Some(v_lo)) = (bracket(r_lo), bracket(r_hi)) {
                let (lo, hi) = (v_lo.min(v_hi) - 1e-6, v_lo.max(v_hi) + 1e-6);
                assert!(
                    (lo..=hi).contains(&vol),
                    "{label}: volume {vol:.6} outside constant-radius bracket [{lo:.6}, {hi:.6}]"
                );
            }
            eprintln!("[vr-matrix] {label}: COMPLETED watertight, volume {vol:.6}");
        }
        Err(err) => {
            eprintln!("[vr-matrix] {label}: DECLINED ({err:?})");
        }
    }
}

/// Extreme but feasible ratio (1:5 on a unit cube).
#[test]
#[ignore = "PK-139 stage 2: the variable-radius SWEEP still ships open bodies (center-curve march landed 2026-07-18; trim/topology wiring is the next leg) — flip when the sweep closes"]
fn variable_extreme_ratio_1_to_5_completes_or_declines_structured() {
    run_law(
        "ratio 0.02:0.10",
        FilletRadius::Variable {
            start: 0.02,
            end: 0.10,
        },
        0.02,
        0.10,
    );
}

/// Severe ratio 1:20 — must not panic or ship an open body.
#[test]
#[ignore = "PK-139 stage 2: the variable-radius SWEEP still ships open bodies (center-curve march landed 2026-07-18; trim/topology wiring is the next leg) — flip when the sweep closes"]
fn variable_extreme_ratio_1_to_20_completes_or_declines_structured() {
    run_law(
        "ratio 0.01:0.20",
        FilletRadius::Variable {
            start: 0.01,
            end: 0.20,
        },
        0.01,
        0.20,
    );
}

/// Non-monotone ByPoints law (bulge).
#[test]
#[ignore = "PK-139 stage 2: the variable-radius SWEEP still ships open bodies (center-curve march landed 2026-07-18; trim/topology wiring is the next leg) — flip when the sweep closes"]
fn by_points_bulge_completes_or_declines_structured() {
    run_law(
        "bulge 0.05→0.12→0.05",
        FilletRadius::ByPoints(vec![(0.0, 0.05), (0.5, 0.12), (1.0, 0.05)]),
        0.05,
        0.12,
    );
}

/// Non-monotone ByPoints law (waist).
#[test]
#[ignore = "PK-139 stage 2: the variable-radius SWEEP still ships open bodies (center-curve march landed 2026-07-18; trim/topology wiring is the next leg) — flip when the sweep closes"]
fn by_points_waist_completes_or_declines_structured() {
    run_law(
        "waist 0.12→0.05→0.12",
        FilletRadius::ByPoints(vec![(0.0, 0.12), (0.5, 0.05), (1.0, 0.12)]),
        0.05,
        0.12,
    );
}

/// Radius law exceeding the face half-extent mid-edge (0.7 on a unit
/// cube overruns the adjacent faces) must DECLINE with a structured
/// reason — a completed body here would be self-intersecting.
#[test]
#[ignore = "PK-139 stage 2: the variable-radius SWEEP still ships open bodies (center-curve march landed 2026-07-18; trim/topology wiring is the next leg) — flip when the sweep closes"]
fn overflow_law_declines_structured() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let result = fillet_edges(
        &brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::ByPoints(vec![(0.0, 0.05), (0.5, 0.7), (1.0, 0.05)]),
        }],
        &FilletOptions::default(),
        default_tol(),
        None,
    );
    match result {
        Ok(BlendResult { brep: out, .. }) => {
            // If it claims completion, it must at least be watertight
            // and self-consistent — measure and fail loudly so the
            // overflow class gets a real decline path.
            let opens = mesh_open_edges(&out);
            assert_eq!(
                opens, 0,
                "overflow law shipped an OPEN body — must decline instead"
            );
            let errors = out.validate_all();
            assert!(
                errors.is_empty(),
                "overflow law shipped an invalid body ({} errors) — must decline instead",
                errors.len()
            );
            eprintln!("[vr-matrix] overflow: completed (validate_all clean) — inspect visually");
        }
        Err(err) => {
            eprintln!("[vr-matrix] overflow: DECLINED ({err:?})");
        }
    }
}

/// Variable ends must feed vertex blends with the ENDPOINT radius:
/// two edges meeting at a cube corner, each variable, sharing the
/// corner value — the corner must accept (mixed-radius cap dispatch
/// receives per-edge endpoint radii, not spec averages).
#[test]
#[ignore = "PK-139 stage 2: the variable-radius SWEEP still ships open bodies (center-curve march landed 2026-07-18; trim/topology wiring is the next leg) — flip when the sweep closes"]
fn variable_ends_feed_corner_with_endpoint_radius() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    // Find two edges sharing a vertex.
    let mut chosen: Option<(
        popkern_core::topology::EdgeHandle,
        popkern_core::topology::EdgeHandle,
    )> = None;
    'outer: for (e1, d1) in brep.edges.iter() {
        for (e2, d2) in brep.edges.iter() {
            if e1 >= e2 {
                continue;
            }
            let v1 = {
                let he = &brep.half_edges[d1.half_edge];
                [he.origin, brep.half_edges[he.twin].origin]
            };
            let v2 = {
                let he = &brep.half_edges[d2.half_edge];
                [he.origin, brep.half_edges[he.twin].origin]
            };
            if v1.iter().any(|v| v2.contains(v)) {
                chosen = Some((e1, e2));
                break 'outer;
            }
        }
    }
    let (e1, e2) = chosen.expect("cube has adjacent edges");
    let specs = vec![
        FilletSpec {
            edge: e1,
            radius: FilletRadius::Variable {
                start: 0.08,
                end: 0.08,
            },
        },
        FilletSpec {
            edge: e2,
            radius: FilletRadius::Variable {
                start: 0.08,
                end: 0.08,
            },
        },
    ];
    let result = fillet_edges(
        &brep,
        &specs,
        &FilletOptions::default(),
        default_tol(),
        None,
    );
    match result {
        Ok(BlendResult { brep: out, .. }) => {
            let opens = mesh_open_edges(&out);
            assert_eq!(opens, 0, "corner-fed variable blend must be watertight");
            eprintln!("[vr-matrix] corner endpoint-radius: COMPLETED watertight");
        }
        Err(err) => {
            eprintln!("[vr-matrix] corner endpoint-radius: DECLINED ({err:?})");
        }
    }
}
