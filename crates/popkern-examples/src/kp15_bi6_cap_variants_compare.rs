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
            let png_path =
                out.join(format!("kp15_bi6_cap_variant_{label}{suffix}.png"));
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
    let tol = Tolerance::default();
    let pyramid = match popkern_core::topology::BRep::n_gonal_pyramid(6, 1.0, 1.2, tol) {
        Ok(p) => p,
        Err(e) => {
            println!("  pyramid_n6 build failed: {e}");
            return;
        }
    };
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0)
        .map(|(h, _)| h);
    let Some(apex) = apex else {
        println!("  pyramid_n6: apex vertex not found");
        return;
    };
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(h, e)| {
            let he = &pyramid.half_edges[e.half_edge];
            let v0 = he.origin;
            let v1 = pyramid.half_edges[he.twin].origin;
            if v0 == apex || v1 == apex { Some(h) } else { None }
        })
        .collect();
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.10),
        })
        .collect();
    let options = FilletOptions::default();
    match fillet_edges(&pyramid, &specs, &options, tol, None) {
        Ok(result) => {
            let step_path = out.join("reference_pyramid_n6_existing_pipeline.step");
            let doc = step_document(&result.brep, "pyramid_n6_reference");
            if let Err(e) = export_step_file(&doc, &step_path, &StepExportOptions::default()) {
                println!("  reference STEP export failed: {e}");
                return;
            }
            println!("  reference: filleted hexagon-6 pyramid (existing pipeline) -> {}",
                     step_path.display());
            let png_path = out.join("reference_pyramid_n6_existing_pipeline.png");
            let camera = popkern_viz::camera::Camera {
                eye: glam::Vec3::new(2.5, 2.5, 1.8),
                target: glam::Vec3::new(0.0, 0.0, 0.6),
                up: glam::Vec3::new(0.0, 0.0, 1.0),
                fov_y: 30.0_f32.to_radians(),
                aspect: 1.0,
                near: 0.01,
                far: 100.0,
            };
            match popkern_viz::headless::render_to_png(
                &result.brep, &camera, 1200, 1200, &png_path, false,
            ) {
                Ok(_) => println!("    reference PNG -> {}", png_path.display()),
                Err(e) => println!("    reference PNG failed: {e}"),
            }
        }
        Err(e) => println!("  pyramid_n6 fillet failed: {e:?}"),
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
        let row0 = [arc_cp(0.0), arc_cp(1.0 / 3.0), arc_cp(2.0 / 3.0), arc_cp(1.0)];
        // Row 1 lifts up slightly toward the apex (positive z direction).
        let lift = 0.03;
        let lift_cp = |t: f64| -> Point3 {
            let x = (1.0 - t) * dx + t * dx_n;
            let y = (1.0 - t) * dy + t * dy_n;
            Point3::new(r * x, r * y, z_rim + lift)
        };
        let row1 = [lift_cp(0.0), lift_cp(1.0 / 3.0), lift_cp(2.0 / 3.0), lift_cp(1.0)];
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
