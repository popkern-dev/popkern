//! KP15 bi-6 polar cap: visual comparison of four construction variants.
//!
//! Produces one STEP file and one PNG per variant on the same hexagon-6
//! input. Variants differ in two axes:
//!
//! 1. **β reparameterization** (KP15 §3.3 eq (17)): whether to transform
//!    the bi-3 tensor border into `q := b ∘ β` before feeding the
//!    constraint matrix.
//! 2. **F_κ minimisation domain**: whether to use the table's default
//!    RREF-determined free column set (lots of degrees of freedom, lower
//!    F_5) or restrict to [KP15 §3.6 step 4]'s "P" set (paper-faithful
//!    but ~3.4× higher F_5).
//!
//! The 2×2 grid:
//!
//! - **A** — no β, full F_5 domain (current `via_table` default)
//! - **B** — no β, KP15 P-only F_5 domain
//! - **C** — with β, full F_5 domain
//! - **D** — with β, KP15 P-only F_5 domain
//!
//! All four variants satisfy eq (2) G¹ and eq (3) G² to machine zero at
//! every cross-sector edge; they differ in the deep-interior CPs the
//! optimisation picks. Use the rendered output to inspect the cap's
//! shape and curvature in a CAD viewer.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::gregory::{
    Bi6Reparams, CirculantSectorRhs, GregoryPatch, Kp15Bi6GeneratingTable, PolarBoundaryInput,
    RectangularGregoryBiNonic, apply_beta_reparameterization, catmull_clark_star_point,
    f_kappa_quadratic_form,
};
use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::geometry::{Point3, SurfaceGeometry, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepProductNode, UnitContext, export_step_file,
};

fn output_dir() -> PathBuf {
    // Match the convention used by `gregory_step_examples`,
    // `boolean_showcase`, etc.: write into the workspace's
    // top-level `examples/output/` tree (which the sandbox serves
    // at `/files/`), one subdirectory per example binary.
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("examples crate has workspace parent")
        .parent()
        .expect("workspace parent has repo parent")
        .join("examples")
        .join("output")
        .join("kp15_bi6_cap_variants_compare")
}

fn main() {
    let out = output_dir();
    std::fs::create_dir_all(&out).expect("create output dir");

    // Reference: real hexagon-6 pyramid filleted via the existing
    // `fillet_edges` pipeline. This is what the kernel currently
    // produces, with the existing apex cap construction (not via-
    // table). Use as the visual ground-truth for "this is how a
    // filleted apex should look."
    export_real_pyramid_reference(&out);

    // Opt-in via-table variants (`--with-via-table`) — they're cap-
    // only loose BReps that don't tessellate cleanly. Focus is the
    // real pyramid renders above.
    if !std::env::args().any(|a| a == "--with-via-table") {
        return;
    }

    // Build the via-table comparison input: a regular hexagon-6 polar
    // config with a symmetric bi-3 tensor border that mimics what a
    // pyramid apex fillet would feed into the cap construction.
    let n = 6;
    let (polar, bi3_per_sector) = hexagon6_polar_and_tensor_border();
    let params = Bi6Reparams::for_valence(n);

    // CirculantSectorRhs from the elevated bi-3 (no β).
    let sectors_no_beta: Vec<CirculantSectorRhs> =
        bi3_per_sector.iter().map(sector_rhs_from_bi3).collect();
    // CirculantSectorRhs from β-reparameterized bi-3.
    let sectors_beta: Vec<CirculantSectorRhs> = bi3_per_sector
        .iter()
        .map(|b| {
            let q = apply_beta_reparameterization(b, &params);
            sector_rhs_from_q_rows(q)
        })
        .collect();

    // Star point (shared across all variants since it depends only on
    // the cap-net, not the tensor border).
    let c7 = polar.central_vertex;
    let mut c5s = Vec::with_capacity(n);
    let mut c6s = Vec::with_capacity(n);
    for spoke in &polar.spokes {
        c5s.push(eval_degree7_bezier(&spoke.boundary_cps, 1.0 / 3.0));
        c6s.push(eval_degree7_bezier(&spoke.boundary_cps, 2.0 / 3.0));
    }
    let star = catmull_clark_star_point(n, c7, &c5s, &c6s);

    let q_local = f_kappa_quadratic_form(5);
    let table_default = Kp15Bi6GeneratingTable::build_with_identifications(n, params);
    let table_kp15_p = Kp15Bi6GeneratingTable::build_with_kp15_p_set_as_free(n, params);

    // Ridge λ on A/C bounds the unconstrained deep-interior corner CPs
    // that the full-F_5 SPD solve would otherwise let balloon. On B/D
    // (KP15 P-only) ridge λ is a no-op because the extras are already
    // pinned at zero — included symmetrically for clarity.
    let lambda_full = 1e-3;
    let lambda_kp15_p = 0.0;

    let (a_grids, a_f) = table_default
        .populate_cap_grids_from_input_with_f_kappa(&sectors_no_beta, star, &q_local, lambda_full)
        .expect("A: no β, full F_5");
    let (b_grids, b_f) = table_kp15_p
        .populate_cap_grids_from_input_with_f_kappa(&sectors_no_beta, star, &q_local, lambda_kp15_p)
        .expect("B: no β, KP15 P-only");
    let (c_grids, c_f) = table_default
        .populate_cap_grids_from_input_with_f_kappa(&sectors_beta, star, &q_local, lambda_full)
        .expect("C: with β, full F_5");
    let (d_grids, d_f) = table_kp15_p
        .populate_cap_grids_from_input_with_f_kappa(&sectors_beta, star, &q_local, lambda_kp15_p)
        .expect("D: with β, KP15 P-only");

    for (label, bi6_grids, fvals) in [
        ("A_noBeta_fullF5", &a_grids, &a_f),
        ("B_noBeta_kp15Ponly", &b_grids, &b_f),
        ("C_withBeta_fullF5", &c_grids, &c_f),
        ("D_withBeta_kp15Ponly", &d_grids, &d_f),
    ] {
        let bi9_grids: Vec<[[Point3; 10]; 10]> = bi6_grids
            .iter()
            .map(|g| elevate_bi6_grid_to_bi9_simple(*g))
            .collect();
        let brep = build_cap_brep(&bi9_grids).expect("build cap BRep");

        let step_path = out.join(format!("kp15_bi6_cap_variant_{label}.step"));
        let doc = step_document(&brep, label);
        export_step_file(&doc, &step_path, &StepExportOptions::default()).expect("STEP export");
        println!(
            "  {label}: F_5 per coord = ({:.3e}, {:.3e}, {:.3e})  →  {}",
            fvals[0],
            fvals[1],
            fvals[2],
            step_path.display()
        );

        let camera = render_camera();
        // Render twice per variant: one shaded (no wireframe) for the
        // cap's actual shape, one with wireframe so the BB control mesh
        // is visible if you want to inspect CP placement.
        for (suffix, wireframe) in [("_shaded", false), ("_wireframe", true)] {
            let png_path = out.join(format!("kp15_bi6_cap_variant_{label}{suffix}.png"));
            match popkern_viz::headless::render_to_png(
                &brep, &camera, 1200, 1200, &png_path, wireframe,
            ) {
                Ok(_) => println!("    PNG{suffix} -> {}", png_path.display()),
                Err(e) => println!("    PNG{suffix} render failed: {e}"),
            }
        }
    }
}

/// Real hexagon-6 pyramid filleted at the apex via the existing
/// fillet pipeline. The output BRep includes the full topology —
/// apex cap, fillet spheres, slope faces, base — and tessellates
/// cleanly. Use this as the visual reference for "what a filleted
/// pyramid actually looks like" so the via-table experiments
/// below can be judged in context.
fn export_real_pyramid_reference(out: &Path) {
    // Just pentagon n=5 sharp h=0.8 r=0.05 — the case that produces a
    // visually recognisable apex cap (RollingBallSphere fallback
    // with G1 ~ 0.86 — fails the 50mrad commercial gate but is the
    // best the kernel currently produces). The hexagon n=6 case and
    // taller pyramids all decline polar Gregory and produce visually
    // worse caps.
    println!("\n=== pentagon n=5, h=0.8, r=0.05 ===");
    run_pyramid_fillet_reference(out, 5, 0.8, "pentagon_n5_sharp");
}

fn run_pyramid_fillet_reference(out: &Path, n: usize, height: f64, label: &str) {
    use popkern_core::topology::{
        BrepExportQualityOptions, BrepTopologyOutputPolicy, BrepTrimRepairOptions,
        analyze_brep_export_quality,
    };
    use popkern_fillet::ApexConstructionDiagnostic;
    let tol = Tolerance::default();
    // Use the same geometry as `kp_polar_commercial_gates.rs::polar_
    // passes_commercial_gates_on_pentagon5_uniform`: base radius 1.0,
    // height 0.8 (sharp apex), fillet radius 0.05.
    let pyramid = match popkern_core::topology::BRep::n_gonal_pyramid(n, 1.0, height, tol) {
        Ok(p) => p,
        Err(e) => {
            println!("  {label} build failed: {e}");
            return;
        }
    };
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0)
        .map(|(h, _)| h);
    let Some(apex) = apex else {
        println!("  {label}: apex vertex not found");
        return;
    };
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(h, e)| {
            let he = &pyramid.half_edges[e.half_edge];
            let v0 = he.origin;
            let v1 = pyramid.half_edges[he.twin].origin;
            if v0 == apex || v1 == apex {
                Some(h)
            } else {
                None
            }
        })
        .collect();
    // Skip the radius sweep — just run the single config.
    if false {
        for radius in &[0.05_f64] {
            let radius = *radius;
            let specs: Vec<FilletSpec> = apex_edges
                .iter()
                .map(|&edge| FilletSpec {
                    edge,
                    radius: FilletRadius::Constant(radius),
                })
                .collect();
            let options = FilletOptions {
                use_rational_gregory: true,
                use_deg5_gregory: true,
                use_biquintic_full_gregory: true,
                use_biseptic_gregory: true,
                use_two_stripe_gregory: true,
                use_kp_binonic_gregory: true,
                use_kp_polar_gregory: true,
                ..FilletOptions::default()
            };
            match fillet_edges(&pyramid, &specs, &options, tol, None) {
                Ok(r) => {
                    let mut worst_g1 = 0.0_f64;
                    let mut decline_codes: Vec<String> = Vec::new();
                    for d in &r.repair_report.apex_diagnostics {
                        match d {
                            ApexConstructionDiagnostic::CapSupportContinuity(c) => {
                                worst_g1 = worst_g1.max(c.worst_g1_coplanarity);
                            }
                            ApexConstructionDiagnostic::ExactBoundaryCap(c)
                                if c.code.contains("decline") =>
                            {
                                decline_codes.push(format!("{}/{}", c.code, c.stage));
                            }
                            _ => {}
                        }
                    }
                    let route = r
                        .repair_report
                        .apex_routes
                        .first()
                        .map(|x| format!("{:?}", x.route))
                        .unwrap_or_else(|| "?".into());
                    let declines = if decline_codes.is_empty() {
                        "no declines".to_string()
                    } else {
                        format!(
                            "{} declines, first = {}",
                            decline_codes.len(),
                            decline_codes[0]
                        )
                    };
                    println!(
                        "    r={radius:.2}: route={route}, worst G1 = {worst_g1:.3e}, {declines}"
                    );
                }
                Err(e) => println!("    r={radius:.2}: FAILED: {e:?}"),
            }
        }
    }

    // Match `polar_passes_commercial_gates_on_pentagon5_uniform`:
    // uniform 0.05 radius across all spokes.
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();
    // Sweep apex routing options to find which produces the best
    // cap. Report worst G1 coplanarity per route; pick the winner.
    let full_gregory_stack = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    let option_variants: Vec<(&str, FilletOptions)> = vec![
        ("default_route", FilletOptions::default()),
        (
            "direct_apex",
            FilletOptions {
                use_direct_apex_construction: true,
                ..FilletOptions::default()
            },
        ),
        ("full_gregory_stack", full_gregory_stack),
    ];
    // Skip the routing sweep — just use the full Gregory stack
    // directly (which is what `kp_polar_commercial_gates` does).
    let (best_name, options) = option_variants
        .iter()
        .find(|(n, _)| *n == "full_gregory_stack")
        .map(|(n, o)| (*n, o.clone()))
        .expect("full_gregory_stack present");
    let result =
        fillet_edges(&pyramid, &specs, &options, tol, None).expect("pentagon n=5 fillet succeeds");
    let mut best_g1 = 0.0_f64;
    for d in &result.repair_report.apex_diagnostics {
        if let ApexConstructionDiagnostic::CapSupportContinuity(c) = d {
            best_g1 = best_g1.max(c.worst_g1_coplanarity);
        }
    }
    println!("  Using {best_name}: worst G1 coplanarity = {best_g1:.3e} (gate = 5.000e-2)");

    println!("  pyramid_n6 fillet returned:");
    println!(
        "    {} corner patches, {} faces, {} edges, {} vertices",
        result.corner_patches.len(),
        result.brep.faces.len(),
        result.brep.edges.len(),
        result.brep.vertices.len()
    );
    // Dump ALL apex diagnostics so we know what actually built the cap.
    println!(
        "  All apex_diagnostics ({}):",
        result.repair_report.apex_diagnostics.len()
    );
    let mut seen_codes: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for d in &result.repair_report.apex_diagnostics {
        match d {
            ApexConstructionDiagnostic::ExactBoundaryCap(c) => {
                let key = format!("{}/{}", c.code, c.stage);
                if seen_codes.insert(key.clone()) {
                    let detail = c
                        .detail_code
                        .as_deref()
                        .map(|d| {
                            d.lines()
                                .next()
                                .unwrap_or("")
                                .chars()
                                .take(120)
                                .collect::<String>()
                        })
                        .unwrap_or_default();
                    println!("    {key}  {detail}");
                }
            }
            ApexConstructionDiagnostic::CapSupportContinuity(c) => {
                println!(
                    "    cap-support continuity: scope={}, worst G1 = {:.3e}, gate = {:.3e}, {} edges measured",
                    c.measurement_scope, c.worst_g1_coplanarity, c.gate, c.edges_measured
                );
            }
            _ => {}
        }
    }
    // List the corner-patch construction method.
    for (i, p) in result.corner_patches.iter().enumerate() {
        println!("    corner patch {i}: method = {:?}", p.method);
    }

    // Apex/cap diagnostics — what the kernel reports about the
    // cap's quality and acceptance.
    let mut apex_route_count = 0usize;
    let mut cap_continuity_count = 0usize;
    let mut worst_g1_co = 0.0_f64;
    let mut declines: Vec<String> = Vec::new();
    for diag in &result.repair_report.apex_diagnostics {
        match diag {
            ApexConstructionDiagnostic::ExactBoundaryCap(c) => {
                if c.code.contains("declined") || c.code.contains("rollback") {
                    declines.push(format!("{}: {}", c.code, c.stage));
                }
            }
            ApexConstructionDiagnostic::CapSupportContinuity(c) => {
                cap_continuity_count += 1;
                if c.worst_g1_coplanarity > worst_g1_co {
                    worst_g1_co = c.worst_g1_coplanarity;
                }
            }
            _ => {}
        }
    }
    for route in &result.repair_report.apex_routes {
        apex_route_count += 1;
        if let Some(reason) = &route.fallback_reason {
            println!("    apex route {:?} fallback: {reason}", route.route);
        } else {
            println!("    apex route {:?} accepted", route.route);
        }
    }
    println!(
        "    {apex_route_count} apex routes attempted, {cap_continuity_count} continuity reports, worst G1 coplanarity = {worst_g1_co:.3e}"
    );
    for d in declines.iter().take(10) {
        println!("    decline: {d}");
    }

    // Topology quality metrics (manifold edges, pcurve deviation,
    // export-readiness).
    let quality = analyze_brep_export_quality(
        &result.brep,
        BrepExportQualityOptions {
            topology_output: BrepTopologyOutputPolicy::OpenSheet,
            trim: BrepTrimRepairOptions {
                endpoint_tolerance: 1e-6,
                sampled_tolerance: 1e-6,
                samples: 9,
                validate_output: true,
            },
            require_edge_curves: true,
            allow_collapsed_trim_edge_curves: true,
            max_pcurve_deviation: Some(1e-6),
            max_edge_surface_deviation: Some(1e-6),
            validate_brep: true,
        },
    );
    println!(
        "    quality: readiness = {:?}, edges = {}, boundary edges = {}",
        quality.readiness, quality.edges, quality.boundary_edges
    );
    println!(
        "    missing edge curves: {}, pcurve blockers: {}, max pcurve dev: {:.3e}",
        quality.missing_edge_curves, quality.pcurve_blockers, quality.max_pcurve_deviation
    );
    println!(
        "    stale pcurves: {}, degenerate topology: {}, collapsed trim missing curves: {}",
        quality.authority.stale_pcurves,
        quality.authority.degenerate_topology,
        quality.collapsed_trim_missing_edge_curves
    );

    let step_path = out.join(format!("reference_{label}.step"));
    let doc = step_document(&result.brep, label);
    if let Err(e) = export_step_file(&doc, &step_path, &StepExportOptions::default()) {
        println!("  {label} STEP export failed: {e}");
        return;
    }
    println!("    STEP -> {}", step_path.display());
    let png_path = out.join(format!("reference_{label}.png"));
    let camera = popkern_viz::camera::Camera {
        eye: glam::Vec3::new(2.0, 2.0, 1.4),
        target: glam::Vec3::new(0.0, 0.0, 0.4),
        up: glam::Vec3::new(0.0, 0.0, 1.0),
        fov_y: 30.0_f32.to_radians(),
        aspect: 1.0,
        near: 0.01,
        far: 100.0,
    };
    match popkern_viz::headless::render_to_png(&result.brep, &camera, 1200, 1200, &png_path, false)
    {
        Ok(_) => println!("    PNG -> {}", png_path.display()),
        Err(e) => println!("    PNG failed: {e}"),
    }

    // Also render close-up on just the apex region.
    let apex_png = out.join(format!("reference_{label}_apex_closeup.png"));
    let apex_camera = popkern_viz::camera::Camera {
        eye: glam::Vec3::new(0.25, 0.25, 1.0),
        target: glam::Vec3::new(0.0, 0.0, 0.78),
        up: glam::Vec3::new(0.0, 0.0, 1.0),
        fov_y: 20.0_f32.to_radians(),
        aspect: 1.0,
        near: 0.01,
        far: 10.0,
    };
    match popkern_viz::headless::render_to_png(
        &result.brep,
        &apex_camera,
        1200,
        1200,
        &apex_png,
        false,
    ) {
        Ok(_) => println!("    apex close-up PNG -> {}", apex_png.display()),
        Err(e) => println!("    apex close-up failed: {e}"),
    }
}

// ─── Input fixtures ─────────────────────────────────────────────────

/// Build a symmetric hexagon-6 polar input that mimics what
/// `BRep::n_gonal_pyramid(6, 1.0, 1.2, ...)` followed by a uniform
/// fillet of radius 0.1 at the apex would produce internally.
///
/// All sectors have identical geometry (no per-sector asymmetry).
/// The tensor border is a small rounded patch near the apex sitting
/// at height z = 1.0 (apex would be at z = 1.2 on the unfilleted
/// pyramid). The "rim" of the cap sits at radius ≈ 0.15 from the
/// apex axis, matching what the rolling-ball fillet would carve.
fn hexagon6_polar_and_tensor_border() -> (PolarBoundaryInput, Vec<[[Point3; 4]; 3]>) {
    use popkern_core::geometry::gregory::{OuterCornerData, SpokeBoundaryData, TangentPlane};
    let n = 6usize;
    let apex = Point3::new(0.0, 0.0, 1.0);
    // The cap's local tangent plane at the apex: roughly horizontal
    // since the cap is a smooth dome.
    let tp_horizontal =
        TangentPlane::from_basis(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap();
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut face_tangent_planes = Vec::with_capacity(n);
    let mut per_spoke_radii = Vec::with_capacity(n);
    let mut bi3_rows = Vec::with_capacity(n);

    // Cap scale: a small region near the apex that the fillet would
    // smooth. Rim radius 0.2 (radial distance from apex axis at the
    // base of the cap), cap depth 0.1 (how far below the apex the
    // rim sits).
    let rim_radius: f64 = 0.2;
    let cap_depth: f64 = 0.1;
    let fillet_radius = 0.1;

    for k in 0..n {
        let theta = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
        let dx = theta.cos();
        let dy = theta.sin();
        // The outer-rim corner where this sector ends, in 3D.
        let rim_point = Point3::new(rim_radius * dx, rim_radius * dy, apex.z - cap_depth);
        // Spoke boundary (degree-7 Bezier from apex outward to rim).
        let mut boundary_cps = [Point3::new(0.0, 0.0, 0.0); 8];
        for (i, cp) in boundary_cps.iter_mut().enumerate() {
            let t = i as f64 / 7.0;
            *cp = Point3::new(
                t * rim_point.x,
                t * rim_point.y,
                apex.z + t * (rim_point.z - apex.z),
            );
        }
        spokes.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c: Vec3::new(-dy, dx, 0.0),
            xtan_b_at_c: Vec3::new(dy, -dx, 0.0),
            xtan_a_at_v: Vec3::new(-dy, dx, 0.0),
            xtan_b_at_v: Vec3::new(dy, -dx, 0.0),
            radius_at_c: fillet_radius,
            radius_at_v: fillet_radius,
            edge_id: format!("edge_{k}"),
        });
        outer_corners.push(OuterCornerData {
            position: rim_point,
            tangent_plane: tp_horizontal,
        });
        face_tangent_planes.push(tp_horizontal);
        per_spoke_radii.push(fillet_radius);

        // Bi-3 tensor border for this sector. Row 0 = position curve
        // along the rim arc from this corner to the next. Row 1 =
        // first cross-derivative (lifting up toward the apex).
        // Row 2 = second cross-derivative (small, for G²).
        // Symmetric across all sectors — no per-sector z_mul.
        let theta_next = (2.0 * std::f64::consts::PI) * ((k + 1) as f64) / (n as f64);
        let dx_n = theta_next.cos();
        let dy_n = theta_next.sin();
        let r = rim_radius;
        let z_rim = apex.z - cap_depth;
        // 4 CPs along the rim arc from sector k's corner to (k+1)'s.
        let arc_cp = |t: f64| -> Point3 {
            let x = (1.0 - t) * dx + t * dx_n;
            let y = (1.0 - t) * dy + t * dy_n;
            Point3::new(r * x, r * y, z_rim)
        };
        let row0 = [
            arc_cp(0.0),
            arc_cp(1.0 / 3.0),
            arc_cp(2.0 / 3.0),
            arc_cp(1.0),
        ];
        // Row 1 lifts up slightly toward the apex (positive z direction).
        let lift = 0.03;
        let lift_cp = |t: f64| -> Point3 {
            let x = (1.0 - t) * dx + t * dx_n;
            let y = (1.0 - t) * dy + t * dy_n;
            Point3::new(r * x, r * y, z_rim + lift)
        };
        let row1 = [
            lift_cp(0.0),
            lift_cp(1.0 / 3.0),
            lift_cp(2.0 / 3.0),
            lift_cp(1.0),
        ];
        // Row 2: small second derivative.
        let lift2 = 0.04;
        let lift2_cp = |t: f64| -> Point3 {
            let x = (1.0 - t) * dx + t * dx_n;
            let y = (1.0 - t) * dy + t * dy_n;
            Point3::new(r * x, r * y, z_rim + lift2)
        };
        let row2 = [
            lift2_cp(0.0),
            lift2_cp(1.0 / 3.0),
            lift2_cp(2.0 / 3.0),
            lift2_cp(1.0),
        ];
        bi3_rows.push([row0, row1, row2]);
    }

    let polar = PolarBoundaryInput {
        central_vertex: apex,
        central_tangent_plane: tp_horizontal,
        spokes,
        outer_corners,
        face_tangent_planes,
        per_spoke_radii,
    };
    (polar, bi3_rows)
}

fn sector_rhs_from_bi3(b: &[[Point3; 4]; 3]) -> CirculantSectorRhs {
    let row0 = elevate_bi3_curve_to_bi6(&b[0]);
    let row1 = elevate_bi3_curve_to_bi6(&b[1]);
    let row2 = elevate_bi3_curve_to_bi6(&b[2]);
    CirculantSectorRhs {
        p_00: row0[0],
        p_10: row0[1],
        p_20: row0[2],
        p_40: row0[4],
        p_50: row0[5],
        p_60: row0[6],
        p_01: row1[0],
        p_11: row1[1],
        p_02: row2[0],
        p_42: row2[4],
        p_52: row2[5],
        p_62: row2[6],
    }
}

fn sector_rhs_from_q_rows(q: [[Point3; 7]; 3]) -> CirculantSectorRhs {
    CirculantSectorRhs {
        p_00: q[0][0],
        p_10: q[0][1],
        p_20: q[0][2],
        p_40: q[0][4],
        p_50: q[0][5],
        p_60: q[0][6],
        p_01: q[1][0],
        p_11: q[1][1],
        p_02: q[2][0],
        p_42: q[2][4],
        p_52: q[2][5],
        p_62: q[2][6],
    }
}

// ─── BRep construction ──────────────────────────────────────────────

/// Build a BRep with one Gregory face per cap sector. Each face has a
/// 4-edge quadrilateral outer loop linking the bi-9 grid's four corner
/// CPs. Topology is per-sector (not shared across sectors); the result
/// is an open shell suitable for STEP export and PNG rendering.
fn build_cap_brep(bi9_grids: &[[[Point3; 10]; 10]]) -> Result<BRep, String> {
    let mut builder = BRepBuilder::new(Tolerance::default());
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    for grid in bi9_grids {
        let v_00 = builder.add_vertex(grid[0][0]);
        let v_90 = builder.add_vertex(grid[9][0]);
        let v_99 = builder.add_vertex(grid[9][9]);
        let v_09 = builder.add_vertex(grid[0][9]);

        let (_, he1, _) = builder.add_edge(v_00, v_90);
        let (_, he2, _) = builder.add_edge(v_90, v_99);
        let (_, he3, _) = builder.add_edge(v_99, v_09);
        let (_, he4, _) = builder.add_edge(v_09, v_00);

        let face = builder.add_face(shell);
        builder
            .add_loop(&[he1, he2, he3, he4], face)
            .map_err(|e| format!("loop: {e:?}"))?;

        let surface = SurfaceOnFace {
            geometry: SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(Box::new(
                RectangularGregoryBiNonic { grid: *grid },
            ))),
            orientation: true,
        };
        builder.set_face_surface(face, surface);
    }

    Ok(builder.build_unchecked())
}

fn step_document(brep: &BRep, label: &str) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: format!("KP15Bi6CapVariant_{label}"),
            id: format!("KP15-CAP-{label}"),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    }
}

// ─── Math helpers ───────────────────────────────────────────────────

fn elevate_bi3_curve_to_bi6(b: &[Point3; 4]) -> [Point3; 7] {
    let mut cur: Vec<Point3> = b.to_vec();
    for _ in 0..3 {
        let d = cur.len() - 1;
        let dp1 = (d + 1) as f64;
        let mut next = Vec::with_capacity(d + 2);
        next.push(cur[0]);
        for i in 1..=d {
            let frac = (i as f64) / dp1;
            let prev = cur[i - 1];
            let this = cur[i];
            next.push(Point3::new(
                frac * prev.x + (1.0 - frac) * this.x,
                frac * prev.y + (1.0 - frac) * this.y,
                frac * prev.z + (1.0 - frac) * this.z,
            ));
        }
        next.push(cur[d]);
        cur = next;
    }
    core::array::from_fn(|i| cur[i])
}

fn elevate_bi6_grid_to_bi9_simple(bi6: [[Point3; 7]; 7]) -> [[Point3; 10]; 10] {
    // First elevate each row from 7 → 10 CPs (degree 6 → 9 in u).
    let mut row_elevated: [[Point3; 10]; 7] = [[Point3::new(0.0, 0.0, 0.0); 10]; 7];
    for (j, row_slot) in row_elevated.iter_mut().enumerate() {
        let row: [Point3; 7] = core::array::from_fn(|i| bi6[i][j]);
        *row_slot = elevate_curve_7_to_10(&row);
    }
    // Then elevate each column from 7 → 10 CPs (degree 6 → 9 in v).
    let mut out: [[Point3; 10]; 10] = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for (i, col_slot) in out.iter_mut().enumerate() {
        let col: [Point3; 7] = core::array::from_fn(|j| row_elevated[j][i]);
        *col_slot = elevate_curve_7_to_10(&col);
    }
    out
}

fn elevate_curve_7_to_10(cps: &[Point3; 7]) -> [Point3; 10] {
    let mut cur: Vec<Point3> = cps.to_vec();
    for _ in 0..3 {
        let d = cur.len() - 1;
        let dp1 = (d + 1) as f64;
        let mut next = Vec::with_capacity(d + 2);
        next.push(cur[0]);
        for i in 1..=d {
            let frac = (i as f64) / dp1;
            let prev = cur[i - 1];
            let this = cur[i];
            next.push(Point3::new(
                frac * prev.x + (1.0 - frac) * this.x,
                frac * prev.y + (1.0 - frac) * this.y,
                frac * prev.z + (1.0 - frac) * this.z,
            ));
        }
        next.push(cur[d]);
        cur = next;
    }
    core::array::from_fn(|i| cur[i])
}

fn eval_degree7_bezier(cps: &[Point3; 8], t: f64) -> Point3 {
    let mut work: [Point3; 8] = *cps;
    for r in 1..=7 {
        for i in 0..(8 - r) {
            let a = work[i];
            let b = work[i + 1];
            work[i] = Point3::new(
                (1.0 - t) * a.x + t * b.x,
                (1.0 - t) * a.y + t * b.y,
                (1.0 - t) * a.z + t * b.z,
            );
        }
    }
    work[0]
}

// ─── Render camera ──────────────────────────────────────────────────

fn render_camera() -> popkern_viz::camera::Camera {
    // Cap sits near the pyramid apex at (0, 0, 1.0) with rim
    // radius 0.2. Position the camera looking down-and-in so the
    // dome is visible from a 3/4 angle, like a typical CAD viewer
    // default for inspecting an apex fillet.
    popkern_viz::camera::Camera {
        eye: glam::Vec3::new(0.6, 0.6, 1.3),
        target: glam::Vec3::new(0.0, 0.0, 0.95),
        up: glam::Vec3::new(0.0, 0.0, 1.0),
        fov_y: 35.0_f32.to_radians(),
        aspect: 1.0,
        near: 0.01,
        far: 10.0,
    }
}
