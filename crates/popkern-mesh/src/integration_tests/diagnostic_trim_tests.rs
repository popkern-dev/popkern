//! Diagnostic data collection for trim regression on showcase examples.

use popkern_core::geometry::Transform3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::TessellationOptions;
use crate::tessellate::tessellate_face;
use crate::trim::{build_trim_loops, trim_loop_area};

/// Return a human-readable name for a [`SurfaceGeometry`] variant.
fn surface_type_name(s: &SurfaceGeometry) -> &'static str {
    match s {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
    }
}

/// Count the number of half-edges in a loop.
fn count_loop_half_edges(brep: &BRep, loop_handle: popkern_core::topology::LoopHandle) -> usize {
    let Some(l) = brep.loops.get(loop_handle) else {
        return 0;
    };
    let start = l.half_edge;
    if !brep.half_edges.contains_key(start) {
        return 0;
    }
    let mut count = 0usize;
    let mut cur = start;
    loop {
        count += 1;
        cur = brep.half_edges[cur].next;
        if cur == start {
            break;
        }
        if count > 10000 {
            break; // safety
        }
    }
    count
}

/// Count how many half-edges in a loop have pcurves.
fn count_pcurves_in_loop(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> (usize, usize) {
    let Some(l) = brep.loops.get(loop_handle) else {
        return (0, 0);
    };
    let start = l.half_edge;
    if !brep.half_edges.contains_key(start) {
        return (0, 0);
    }
    let mut with = 0usize;
    let mut without = 0usize;
    let mut cur = start;
    loop {
        if brep.half_edges[cur].pcurve.is_some() {
            with += 1;
        } else {
            without += 1;
        }
        cur = brep.half_edges[cur].next;
        if cur == start {
            break;
        }
        if with + without > 10000 {
            break;
        }
    }
    (with, without)
}

/// Diagnose a single BRep model: prints per-face trim and tessellation data.
fn diagnose_model(name: &str, brep: &BRep) {
    println!("\n============================================================");
    println!("  MODEL: {name}");
    println!(
        "  Faces: {}, Edges: {}, Vertices: {}",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
    println!("============================================================");

    let opts = TessellationOptions::default();

    for (fi, (fh, face)) in brep.faces.iter().enumerate() {
        let surface_on_face = match &face.surface {
            Some(s) => s,
            None => {
                println!("  Face {fi}: NO SURFACE GEOMETRY -- skipping");
                continue;
            }
        };
        let surface = &surface_on_face.geometry;
        let stype = surface_type_name(surface);

        let outer_he_count = count_loop_half_edges(brep, face.outer_loop);
        let inner_loop_count = face.inner_loops.len();

        println!(
            "\n  Face {fi}: surface={stype}, outer_loop_he={outer_he_count}, \
                  inner_loops={inner_loop_count}"
        );

        // Pcurve counts for outer loop.
        let (pc_with, pc_without) = count_pcurves_in_loop(brep, face.outer_loop);
        println!("    Outer loop pcurves: {pc_with} with, {pc_without} without");

        // Inner loop pcurve counts.
        for (il_idx, &il) in face.inner_loops.iter().enumerate() {
            let il_he = count_loop_half_edges(brep, il);
            let (il_pc_w, il_pc_wo) = count_pcurves_in_loop(brep, il);
            println!(
                "    Inner loop {il_idx}: he={il_he}, pcurves: \
                      {il_pc_w} with, {il_pc_wo} without"
            );
        }

        // Build trim loops.
        let trim_loops = build_trim_loops(brep, fh, surface, 0.01);
        let area = trim_loop_area(&trim_loops);
        println!(
            "    Trim loops: count={}, outer_area={area:.6}",
            trim_loops.len()
        );

        for (tl_idx, tl) in trim_loops.iter().enumerate() {
            let tl_area = crate::trim::shoelace_area(&tl.points);
            println!(
                "      loop {tl_idx}: is_outer={}, points={}, area={tl_area:.6}",
                tl.is_outer,
                tl.points.len()
            );
        }

        // Tessellate this face.
        match tessellate_face(brep, fh, &opts) {
            Ok(fm) => {
                let nv = fm.vertices.len();
                let nt = fm.triangles.len();

                // Compute bounding box of vertex positions.
                let mut x_min = f64::INFINITY;
                let mut x_max = f64::NEG_INFINITY;
                let mut y_min = f64::INFINITY;
                let mut y_max = f64::NEG_INFINITY;
                let mut z_min = f64::INFINITY;
                let mut z_max = f64::NEG_INFINITY;

                for v in &fm.vertices {
                    x_min = x_min.min(v.position.x);
                    x_max = x_max.max(v.position.x);
                    y_min = y_min.min(v.position.y);
                    y_max = y_max.max(v.position.y);
                    z_min = z_min.min(v.position.z);
                    z_max = z_max.max(v.position.z);
                }

                println!("    Tessellation: verts={nv}, tris={nt}");
                println!(
                    "    Bounds: x=[{x_min:.4}, {x_max:.4}] \
                          y=[{y_min:.4}, {y_max:.4}] z=[{z_min:.4}, {z_max:.4}]"
                );

                // Flag extreme vertices.
                let mut extreme_count = 0usize;
                for (vi, v) in fm.vertices.iter().enumerate() {
                    if v.position.x.abs() > 10.0
                        || v.position.y.abs() > 10.0
                        || v.position.z.abs() > 10.0
                        || v.position.x.is_nan()
                        || v.position.y.is_nan()
                        || v.position.z.is_nan()
                    {
                        if extreme_count < 10 {
                            println!(
                                "    *** EXTREME vertex {vi}: pos=({:.4}, {:.4}, {:.4}) \
                                 uv=({:.6}, {:.6})",
                                v.position.x, v.position.y, v.position.z, v.uv.0, v.uv.1
                            );
                        }
                        extreme_count += 1;
                    }
                }
                if extreme_count > 10 {
                    println!("    ... and {} more extreme vertices", extreme_count - 10);
                }
                if extreme_count > 0 {
                    println!("    TOTAL EXTREME VERTICES: {extreme_count}/{nv}");
                }
            }
            Err(e) => {
                println!("    Tessellation FAILED: {e}");
            }
        }
    }
}

// ── Build each showcase model and diagnose ──────────────────────────

fn build_shelled_box() -> Option<BRep> {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).ok()?;
    let face_to_remove = cube.faces.keys().next()?;
    match popkern_shell::shell(
        &cube,
        0.1,
        popkern_shell::ShellDirection::Inward,
        &[face_to_remove],
        &[],
        &popkern_shell::ShellOptions::default(),
        tol,
    ) {
        Ok(result) => Some(result.brep),
        Err(e) => {
            println!("  Shell failed: {e} -- diagnosing plain cube instead");
            Some(cube)
        }
    }
}

fn build_fillet_showcase() -> Option<BRep> {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).ok()?;
    let edge_handles: Vec<_> = cube.edges.keys().collect();
    let specs: Vec<popkern_fillet::FilletSpec> = edge_handles
        .iter()
        .map(|&edge| popkern_fillet::FilletSpec {
            edge,
            radius: popkern_fillet::FilletRadius::Constant(0.1),
        })
        .collect();
    match popkern_fillet::fillet_edges(
        &cube,
        &specs,
        &popkern_fillet::FilletOptions::default(),
        tol,
        None,
    ) {
        Ok(result) => Some(result.brep),
        Err(e) => {
            println!("  Full fillet failed: {e} -- trying single edge");
            let single = vec![popkern_fillet::FilletSpec {
                edge: edge_handles[0],
                radius: popkern_fillet::FilletRadius::Constant(0.1),
            }];
            match popkern_fillet::fillet_edges(
                &cube,
                &single,
                &popkern_fillet::FilletOptions::default(),
                tol,
                None,
            ) {
                Ok(r) => Some(r.brep),
                Err(e2) => {
                    println!("  Single fillet also failed: {e2} -- diagnosing cube");
                    Some(cube)
                }
            }
        }
    }
}

fn build_boolean_showcase() -> Option<BRep> {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).ok()?;
    let hole = BRep::cylinder(0.2, 2.0, tol).ok()?;
    let hole = hole
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .ok()?;
    match popkern_bool::boolean(
        popkern_bool::BoolOperation::Subtraction,
        cube.clone(),
        hole,
        tol,
    ) {
        Ok(brep) => Some(brep),
        Err(e) => {
            println!("  Boolean subtraction failed: {e} -- diagnosing cube");
            Some(cube)
        }
    }
}

#[test]
fn diagnostic_trim_regression() {
    println!("\n\n========== TRIM REGRESSION DIAGNOSTICS ==========\n");

    // 1. Shelled box
    println!("\n--- Building shelled_box ---");
    if let Some(brep) = build_shelled_box() {
        diagnose_model("shelled_box", &brep);
    } else {
        println!("  FAILED to build shelled_box model");
    }

    // 2. Fillet showcase
    println!("\n--- Building fillet_showcase ---");
    if let Some(brep) = build_fillet_showcase() {
        diagnose_model("fillet_showcase", &brep);
    } else {
        println!("  FAILED to build fillet_showcase model");
    }

    // 3. Boolean showcase
    println!("\n--- Building boolean_showcase ---");
    if let Some(brep) = build_boolean_showcase() {
        // Extra: find face(s) with inner loops and report trim usage.
        let mut found_inner = false;
        for (fh, face) in &brep.faces {
            if !face.inner_loops.is_empty() {
                found_inner = true;
                if let Some(sof) = &face.surface {
                    let surface = &sof.geometry;
                    let trim_loops = build_trim_loops(&brep, fh, surface, 0.01);
                    let inner_count = trim_loops.iter().filter(|tl| !tl.is_outer).count();
                    println!(
                        "  Boolean inner-loop face: total trim loops={}, \
                         inner trim loops={inner_count}",
                        trim_loops.len()
                    );
                    for (idx, tl) in trim_loops.iter().enumerate() {
                        let a = crate::trim::shoelace_area(&tl.points);
                        println!(
                            "    trim loop {idx}: is_outer={}, pts={}, area={a:.6}",
                            tl.is_outer,
                            tl.points.len()
                        );
                    }
                }
            }
        }
        if !found_inner {
            println!("  Boolean result has NO faces with inner loops");
        }

        diagnose_model("boolean_showcase", &brep);
    } else {
        println!("  FAILED to build boolean_showcase model");
    }

    println!("\n========== END DIAGNOSTICS ==========\n");
}

#[test]
fn print_single_fillet_vertices() {
    use popkern_core::geometry::Transform3;
    let tol = popkern_core::tolerance::Tolerance::default();
    let cube = popkern_core::topology::BRep::unit_cube(tol).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let result = popkern_fillet::fillet_edge(&cube, edge, 0.1, tol, None).unwrap();
    eprintln!(
        "Single-edge fillet: {} vertices",
        result.brep.vertices.len()
    );
    for (vh, v) in result.brep.vertices.iter() {
        let p = v.position;
        // Only print non-cube-corner vertices
        let is_corner =
            [0.0, 1.0].contains(&p.x) && [0.0, 1.0].contains(&p.y) && [0.0, 1.0].contains(&p.z);
        if !is_corner {
            eprintln!("  NEW {:?}: ({:.6}, {:.6}, {:.6})", vh, p.x, p.y, p.z);
        }
    }
}

#[test]
fn print_fillet_footprint_positions() {
    use popkern_core::geometry::surface::SurfaceEval;
    let tol = popkern_core::tolerance::Tolerance::default();
    let cube = popkern_core::topology::BRep::unit_cube(tol).unwrap();
    let edge = cube.edges.keys().next().unwrap();

    // Print the edge being filleted
    let he = cube.edges[edge].half_edge;
    let v_start = cube.vertices[cube.half_edges[he].origin].position;
    let twin = cube.half_edges[he].twin;
    let v_end = cube.vertices[cube.half_edges[twin].origin].position;
    eprintln!(
        "Filleting edge from ({:.3},{:.3},{:.3}) to ({:.3},{:.3},{:.3})",
        v_start.x, v_start.y, v_start.z, v_end.x, v_end.y, v_end.z
    );

    let result = popkern_fillet::fillet_edge(&cube, edge, 0.1, tol, None).unwrap();
    // Print all non-corner vertices
    for (_, v) in result.brep.vertices.iter() {
        let p = v.position;
        let is_corner = [0.0_f64, 1.0].iter().any(|&c| (p.x - c).abs() < 1e-6)
            && [0.0_f64, 1.0].iter().any(|&c| (p.y - c).abs() < 1e-6)
            && [0.0_f64, 1.0].iter().any(|&c| (p.z - c).abs() < 1e-6);
        if !is_corner {
            eprintln!("  fillet vtx: ({:.6},{:.6},{:.6})", p.x, p.y, p.z);
        }
    }
}

#[test]
fn trace_fillet_center_curve() {
    use popkern_core::geometry::curve::CurveEval;
    let tol = popkern_core::tolerance::Tolerance::default();
    let cube = popkern_core::topology::BRep::unit_cube(tol).unwrap();
    let edge = cube.edges.keys().next().unwrap();

    // Get the analysis for this edge
    let analyses = popkern_fillet::analysis::analyze_edges(&cube, &[edge], false).unwrap();
    let analysis = &analyses[0];

    eprintln!(
        "Surface 1: {:?}",
        std::mem::discriminant(&analysis.surface_1)
    );
    eprintln!(
        "Surface 2: {:?}",
        std::mem::discriminant(&analysis.surface_2)
    );

    // Compute the center curve
    let ((u1_lo, u1_hi), (v1_lo, v1_hi)) =
        popkern_core::geometry::uv_projection::bounded_parameter_range(
            &cube.face_outer_loop_positions(analysis.face_1),
            &analysis.surface_1,
        );
    let ((u2_lo, u2_hi), (v2_lo, v2_hi)) =
        popkern_core::geometry::uv_projection::bounded_parameter_range(
            &cube.face_outer_loop_positions(analysis.face_2),
            &analysis.surface_2,
        );

    eprintln!(
        "Face 1 UV range: u=[{:.3},{:.3}], v=[{:.3},{:.3}]",
        u1_lo, u1_hi, v1_lo, v1_hi
    );
    eprintln!(
        "Face 2 UV range: u=[{:.3},{:.3}], v=[{:.3},{:.3}]",
        u2_lo, u2_hi, v2_lo, v2_hi
    );

    let edge_mid = {
        let he = &cube.half_edges[cube.edges[edge].half_edge];
        let p1 = cube.vertices[he.origin].position;
        let twin = &cube.half_edges[he.twin];
        let p2 = cube.vertices[twin.origin].position;
        popkern_core::geometry::Point3::new(
            (p1.x + p2.x) * 0.5,
            (p1.y + p2.y) * 0.5,
            (p1.z + p2.z) * 0.5,
        )
    };
    eprintln!(
        "Edge midpoint: ({:.3},{:.3},{:.3})",
        edge_mid.x, edge_mid.y, edge_mid.z
    );

    let cc = popkern_fillet::center_curve::constant_radius_center_curve(
        &analysis.surface_1,
        (u1_lo, u1_hi),
        (v1_lo, v1_hi),
        &analysis.surface_2,
        (u2_lo, u2_hi),
        (v2_lo, v2_hi),
        0.1,
        &edge_mid,
        &tol,
        analysis.orientation_1,
        analysis.orientation_2,
    );

    match cc {
        Ok(result) => {
            let (t0, t1) = result.center_curve.parameter_range();
            for i in 0..5 {
                let t = t0 + (t1 - t0) * (i as f64 / 4.0);
                let p = result.center_curve.point_at(t);
                eprintln!("  center[{i}]: ({:.6},{:.6},{:.6})", p.x, p.y, p.z);
            }
            eprintln!("  footprint_1: {} points", result.footprint_1.len());
            for (i, fp) in result.footprint_1.iter().enumerate().take(5) {
                eprintln!("    fp1[{i}]: ({:.6},{:.6},{:.6})", fp.x, fp.y, fp.z);
            }
            eprintln!("  footprint_2: {} points", result.footprint_2.len());
            for (i, fp) in result.footprint_2.iter().enumerate().take(5) {
                eprintln!("    fp2[{i}]: ({:.6},{:.6},{:.6})", fp.x, fp.y, fp.z);
            }
        }
        Err(e) => eprintln!("  center curve failed: {e}"),
    }
}
