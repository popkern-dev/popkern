//! Step 1 of the post-D1.1 commercial-grade investigation: verify
//! that our cap surface's boundary along the cap/support seam is a
//! degree-6 polynomial in the bi-9 Gregory representation, matching
//! the KP20 paper's prescription (the support cells are bi-6 in
//! `Kp20Bi6Patch`, elevated to bi-9 for storage in
//! `GregoryPatch::RectangularBiNonic`).
//!
//! Approach: for every `RectangularBiNonic` cap/support face in the
//! hexagon6_variable BRep, extract each side's 10 boundary CPs
//! (`grid[0][..]`, `grid[9][..]`, `grid[..][0]`, `grid[..][9]`) and
//! convert from Bernstein basis to monomial basis. A polynomial of
//! true degree `d` has zero monomial coefficients at indices `d+1`
//! through `9`. If our seam boundaries are genuinely degree-6 (as
//! KP20 prescribes), the `t^7`, `t^8`, `t^9` coefficient vectors
//! should be at machine precision.
//!
//! Tagged `#[ignore]` by default — observational probe.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence};

fn run_n_pyramid(n: usize, radii: &[f64]) -> BRep {
    let pyramid = BRep::n_gonal_pyramid(n, 1.0, 0.8, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyramid, n).expect("apex must exist");
    let edges = incident_edges(&pyramid, apex);
    let specs: Vec<_> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let result = fillet_edges(
        &pyramid,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("fillet should not panic");
    result.brep
}

fn run_cube(radii: &[f64]) -> BRep {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let specs: Vec<_> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let result = fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("cube fillet should not panic");
    result.brep
}

fn run_hexagon6_variable() -> BRep {
    run_n_pyramid(6, &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02])
}

/// Binomial coefficient `C(n, k)` for small `n` (≤ 9).
fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let mut num = 1.0_f64;
    let mut den = 1.0_f64;
    for i in 0..k {
        num *= (n - i) as f64;
        den *= (i + 1) as f64;
    }
    num / den
}

/// Convert `n+1` Bezier control points to the monomial coefficient
/// at index `k`. Formula:
///
/// `a_k = Σ_{i=0..=k} (-1)^{k-i} · C(n, i) · C(n-i, k-i) · P_i`
///
/// where `B(t) = Σ_i P_i · C(n,i) · t^i · (1-t)^{n-i}
///             = Σ_k a_k · t^k`.
fn bernstein_to_monomial(cps: &[Point3], k: usize) -> Point3 {
    let n = cps.len() - 1;
    let mut sum = Point3::new(0.0, 0.0, 0.0);
    for (i, cp) in cps.iter().enumerate().take(k.min(n) + 1) {
        let coeff = binomial(n, i) * binomial(n - i, k - i);
        let sign = if (k - i).is_multiple_of(2) { 1.0 } else { -1.0 };
        sum = Point3::new(
            sum.x + sign * coeff * cp.x,
            sum.y + sign * coeff * cp.y,
            sum.z + sign * coeff * cp.z,
        );
    }
    sum
}

fn point_norm(p: Point3) -> f64 {
    (p.x * p.x + p.y * p.y + p.z * p.z).sqrt()
}

fn extract_boundary_row(grid: &[[Point3; 10]; 10], side: BoundarySide) -> [Point3; 10] {
    match side {
        BoundarySide::V0 => grid[0],
        BoundarySide::V1 => grid[9],
        BoundarySide::U0 => core::array::from_fn(|i| grid[i][0]),
        BoundarySide::U1 => core::array::from_fn(|i| grid[i][9]),
    }
}

#[derive(Copy, Clone, Debug)]
enum BoundarySide {
    V0,
    V1,
    U0,
    U1,
}

impl BoundarySide {
    fn label(&self) -> &'static str {
        match self {
            BoundarySide::V0 => "v=0",
            BoundarySide::V1 => "v=1",
            BoundarySide::U0 => "u=0",
            BoundarySide::U1 => "u=1",
        }
    }
}

fn surface_kind_label(geometry: &SurfaceGeometry) -> &'static str {
    match geometry {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(g) => match g {
            GregoryPatch::RectangularBiNonic(_) => "Gregory:RectangularBiNonic",
            GregoryPatch::RectangularBiNonicDualTwist(_) => "Gregory:RectangularBiNonicDualTwist",
            _ => "Gregory:other",
        },
        SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
    }
}

fn report_surface_inventory(label: &str, brep: &BRep) {
    let mut surface_kind_counts: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();
    let mut surfaceless_count = 0usize;
    for (_, face) in brep.faces.iter() {
        match &face.surface {
            Some(s) => {
                *surface_kind_counts
                    .entry(surface_kind_label(&s.geometry))
                    .or_insert(0) += 1
            }
            None => surfaceless_count += 1,
        }
    }
    eprintln!(
        "[{}] total {} faces, {} surfaceless",
        label,
        brep.faces.len(),
        surfaceless_count
    );
    for (kind, count) in &surface_kind_counts {
        eprintln!("    {}: {}", kind, count);
    }
}

#[test]
#[ignore = "fixture-matrix surface inventory probe"]
fn fixture_matrix_surface_kind_inventory() {
    eprintln!("\n## Final-BRep surface inventory across the KP20 polar fixture matrix\n");
    report_surface_inventory("cube_uniform", &run_cube(&[0.10, 0.10, 0.10]));
    report_surface_inventory("cube_mixed", &run_cube(&[0.08, 0.10, 0.12]));
    report_surface_inventory("cube_extreme", &run_cube(&[0.05, 0.10, 0.15]));
    report_surface_inventory(
        "pyramid4_mixed",
        &run_n_pyramid(4, &[0.04, 0.06, 0.05, 0.07]),
    );
    report_surface_inventory(
        "pyramid5_mixed",
        &run_n_pyramid(5, &[0.04, 0.06, 0.05, 0.07, 0.05]),
    );
    report_surface_inventory("hexagon6_uniform", &run_n_pyramid(6, &[0.04; 6]));
    report_surface_inventory(
        "hexagon6_mixed",
        &run_n_pyramid(6, &[0.04, 0.06, 0.05, 0.07, 0.04, 0.06]),
    );
    report_surface_inventory("hexagon6_variable", &run_hexagon6_variable());
}

#[test]
#[ignore = "rank-test probe — run with --ignored --nocapture"]
fn hexagon6_variable_cap_surface_boundary_degree_rank() {
    let brep = run_hexagon6_variable();
    let mut binonic_faces = Vec::new();
    let mut surface_kind_counts: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();
    let mut surfaceless_count = 0usize;
    for (handle, face) in brep.faces.iter() {
        match &face.surface {
            Some(s) => {
                let label = surface_kind_label(&s.geometry);
                *surface_kind_counts.entry(label).or_insert(0) += 1;
                if let SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(b)) = &s.geometry {
                    binonic_faces.push((handle, b.grid));
                }
            }
            None => surfaceless_count += 1,
        }
    }

    eprintln!("\n## hexagon6_variable RectangularBiNonic cap/support face boundary degree rank\n");
    eprintln!("Total faces: {}", brep.faces.len());
    eprintln!("Surfaceless faces: {}", surfaceless_count);
    eprintln!("Surface kind inventory:");
    for (kind, count) in &surface_kind_counts {
        eprintln!("  {}: {}", kind, count);
    }
    eprintln!();
    eprintln!(
        "{} faces with `RectangularBiNonic` surface",
        binonic_faces.len()
    );
    eprintln!();
    eprintln!("| face | side | extent (m) | ‖a_7‖ | ‖a_8‖ | ‖a_9‖ | ratio_7 | ratio_8 | ratio_9 |");
    eprintln!("|---|---|---|---|---|---|---|---|---|");

    let mut worst_ratio_7 = 0.0_f64;
    let mut worst_ratio_8 = 0.0_f64;
    let mut worst_ratio_9 = 0.0_f64;

    for (face, grid) in &binonic_faces {
        for side in [
            BoundarySide::V0,
            BoundarySide::V1,
            BoundarySide::U0,
            BoundarySide::U1,
        ] {
            let cps = extract_boundary_row(grid, side);
            let extent = point_norm(Point3::new(
                cps[9].x - cps[0].x,
                cps[9].y - cps[0].y,
                cps[9].z - cps[0].z,
            ));
            let a7 = point_norm(bernstein_to_monomial(&cps, 7));
            let a8 = point_norm(bernstein_to_monomial(&cps, 8));
            let a9 = point_norm(bernstein_to_monomial(&cps, 9));
            let r7 = if extent > 1e-15 { a7 / extent } else { 0.0 };
            let r8 = if extent > 1e-15 { a8 / extent } else { 0.0 };
            let r9 = if extent > 1e-15 { a9 / extent } else { 0.0 };
            worst_ratio_7 = worst_ratio_7.max(r7);
            worst_ratio_8 = worst_ratio_8.max(r8);
            worst_ratio_9 = worst_ratio_9.max(r9);
            eprintln!(
                "| {:?} | {} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | {:.3e} | {:.3e} |",
                face,
                side.label(),
                extent,
                a7,
                a8,
                a9,
                r7,
                r8,
                r9
            );
        }
    }

    eprintln!();
    eprintln!(
        "Worst ratios (||a_k|| / boundary_extent):  a_7={:.3e}  a_8={:.3e}  a_9={:.3e}",
        worst_ratio_7, worst_ratio_8, worst_ratio_9
    );
    eprintln!();
    eprintln!("Interpretation:");
    eprintln!(
        "  - All three ratios at ~1e-14 to ~1e-12 → boundary is degree-6 in bi-9 form (matches KP20)."
    );
    eprintln!(
        "  - Any ratio > ~1e-6 → boundary carries genuine degree-7+ content; cap construction deviates from KP20."
    );
    eprintln!(
        "  - Mixed signal (some at 1e-14, others larger) → some faces match KP20, others use a non-paper path."
    );
}
