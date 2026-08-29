//! Shared helpers for popkern showcase examples.

use std::path::Path;
use std::sync::Arc;

use popkern_core::topology::BRep;

/// Print a one-shot health card: topology counts, geometry stats,
/// defect categories. Complements `report_non_manifold`. Set
/// `POPKERN_DIAG_VERBOSE=1` for the full per-defect listing.
#[allow(dead_code)]
pub fn report_diagnostics(brep: &BRep) {
    let report = popkern_core::diagnostics::Diagnostics::run(brep);
    print!("{report}");
}

/// For every edge whose tessellation divergence exceeds the env-
/// var-controlled threshold (default 1.0), print a 9-point sample
/// comparison: 3-D curve point, surface_a evaluation via pcurve_a,
/// surface_b evaluation via pcurve_b. Lets us see exactly where
/// things diverge along the curve and on which side the geometry
/// is broken.
#[allow(dead_code)]
pub fn debug_extreme_non_manifold(brep: &BRep) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let threshold = 5.0 * brep.tolerance.linear().max(1e-9);
    let reports = brep.find_tessellation_non_manifold_edges(16, threshold);
    let cutoff = 1.0;
    let extreme: Vec<_> = reports
        .iter()
        .filter(|r| r.max_divergence > cutoff)
        .collect();
    if extreme.is_empty() {
        return;
    }
    println!("  -- extreme non-manifold edges (div > {cutoff}) --");
    for r in &extreme {
        let Some(edge) = brep.edges.get(r.edge) else {
            continue;
        };
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let twin = he.twin;
        let Some(twin_he) = brep.half_edges.get(twin) else {
            continue;
        };
        let pa = he.pcurve.as_ref();
        let pb = twin_he.pcurve.as_ref();
        let face_a = brep.loops[he.loop_].face;
        let face_b = brep.loops[twin_he.loop_].face;
        let surf_a = brep.faces.get(face_a).and_then(|f| f.surface.as_ref());
        let surf_b = brep.faces.get(face_b).and_then(|f| f.surface.as_ref());
        let (Some(pa), Some(pb), Some(sa), Some(sb)) = (pa, pb, surf_a, surf_b) else {
            continue;
        };
        println!(
            "  Edge {:?} (max div={:.3} at param-of-curve)",
            r.edge, r.max_divergence
        );
        let (t0, t1) = (coe.t_start, coe.t_end);
        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        for i in 0..9 {
            let s = i as f64 / 8.0;
            let t = t0 + s * (t1 - t0);
            let p_3d = coe.geometry.point_at(t);
            let ta = pa_lo + s * (pa_hi - pa_lo);
            let tb = pb_lo + s * (pb_hi - pb_lo);
            let uva = pa.point_at(ta);
            let uvb = pb.point_at(tb);
            let p_a = sa.geometry.point_at(uva.x, uva.y);
            let p_b = sb.geometry.point_at(uvb.x, uvb.y);
            println!(
                "    s={:.2}: 3D=({:+.3},{:+.3},{:+.3})  A=({:+.3},{:+.3},{:+.3}) uv=({:+.3},{:+.3})  B=({:+.3},{:+.3},{:+.3}) uv=({:+.3},{:+.3})",
                s,
                p_3d.x,
                p_3d.y,
                p_3d.z,
                p_a.x,
                p_a.y,
                p_a.z,
                uva.x,
                uva.y,
                p_b.x,
                p_b.y,
                p_b.z,
                uvb.x,
                uvb.y,
            );
        }
    }
}

/// Print a count + breakdown of edges where surface_a(pcurve_a) and
/// surface_b(pcurve_b) disagree at interior points — i.e. the edges
/// that OCCT-strict slicers (Bambu Studio, FreeCAD) flag as
/// non-manifold when they tessellate the imported STEP.
///
/// Enable the `non-manifold-verbose` trace topic to also dump every
/// offending edge
/// (handle, surface-pair kind, both faces, location, divergence).
#[allow(dead_code)]
pub fn report_non_manifold(brep: &BRep) {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let threshold = 5.0 * brep.tolerance.linear().max(1e-9);
    let reports = brep.find_tessellation_non_manifold_edges(16, threshold);
    if reports.is_empty() {
        println!("  Non-manifold (tess.): 0");
        return;
    }
    let surf_kind = |sg: &SurfaceGeometry| -> &'static str {
        match sg {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cyl",
            SurfaceGeometry::Sphere { .. } => "Sph",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
            SurfaceGeometry::Gregory(_) => "Gregory",
            SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
        }
    };
    let pair_label = |r: &popkern_core::validation::NonManifoldEdgeReport| -> String {
        let a = brep
            .faces
            .get(r.face_a)
            .and_then(|f| f.surface.as_ref())
            .map(|s| surf_kind(&s.geometry))
            .unwrap_or("?");
        let b = brep
            .faces
            .get(r.face_b)
            .and_then(|f| f.surface.as_ref())
            .map(|s| surf_kind(&s.geometry))
            .unwrap_or("?");
        // Symmetrise (Plane×Cyl == Cyl×Plane)
        if a <= b {
            format!("{a}×{b}")
        } else {
            format!("{b}×{a}")
        }
    };

    // Magnitude buckets.
    let mut buckets: [usize; 4] = [0; 4];
    for r in &reports {
        let div = r.max_divergence;
        if div > 1e-1 {
            buckets[3] += 1;
        } else if div > 1e-2 {
            buckets[2] += 1;
        } else if div > 1e-3 {
            buckets[1] += 1;
        } else {
            buckets[0] += 1;
        }
    }
    println!(
        "  Non-manifold (tess.): {} edges (>1e-1: {}, >1e-2: {}, >1e-3: {}, ≤1e-3: {})",
        reports.len(),
        buckets[3],
        buckets[2],
        buckets[1],
        buckets[0],
    );

    // Surface-pair breakdown — the most actionable view: tells you
    // *which* class of edges is misbehaving (Plane×Cyl drift vs
    // Nurbs×Cyl cap, etc.).
    let mut by_pair: std::collections::BTreeMap<String, (usize, f64)> =
        std::collections::BTreeMap::new();
    for r in &reports {
        let key = pair_label(r);
        let entry = by_pair.entry(key).or_insert((0, 0.0));
        entry.0 += 1;
        if r.max_divergence > entry.1 {
            entry.1 = r.max_divergence;
        }
    }
    let mut pairs: Vec<_> = by_pair.into_iter().collect();
    pairs.sort_by_key(|p| std::cmp::Reverse(p.1.0));
    print!("    by pair: ");
    let mut first = true;
    for (k, (n, max_div)) in &pairs {
        if !first {
            print!(", ");
        }
        first = false;
        print!("{k}={n} (max {:.3})", max_div);
    }
    println!();

    // Top offenders — handle + faces + location for hand-investigation.
    let n_show = reports.len().min(3);
    for r in &reports[..n_show] {
        println!(
            "    worst: {:?} ({}) {:?}|{:?} at ({:+.3},{:+.3},{:+.3}) div={:.4}",
            r.edge,
            pair_label(r),
            r.face_a,
            r.face_b,
            r.location.x,
            r.location.y,
            r.location.z,
            r.max_divergence
        );
    }

    if popkern_core::trace::trace_enabled("non-manifold-verbose") {
        println!("    -- all non-manifold edges --");
        for r in &reports {
            println!(
                "      {:?} ({}) {:?}|{:?} at ({:+.3},{:+.3},{:+.3}) div={:.5}",
                r.edge,
                pair_label(r),
                r.face_a,
                r.face_b,
                r.location.x,
                r.location.y,
                r.location.z,
                r.max_divergence
            );
        }
    }
}
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepProductNode, UnitContext, export_step_file,
};

/// Export a BRep solid as a STEP file with the given product name.
#[allow(dead_code)]
pub fn export(brep: &BRep, name: &str, filename: &str) {
    let output_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("examples")
        .join("output");
    std::fs::create_dir_all(&output_dir).expect("create output dir");
    let path = output_dir.join(filename);

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: name.into(),
            id: format!("{}-001", name.to_uppercase()),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    match export_step_file(&doc, &path, &StepExportOptions::default()) {
        Ok(()) => {
            println!("  Exported: {}", path.display());
            // Validate by re-importing the STEP file.
            validate_step_file(&path);
        }
        Err(e) => eprintln!("  STEP export failed: {e}"),
    }
}

/// Validate a STEP file by re-importing it and checking for errors.
#[allow(dead_code)]
fn validate_step_file(path: &Path) {
    use popkern_io::{StepImportOptions, import_step};

    match std::fs::read(path) {
        Ok(data) => match import_step(&data, &StepImportOptions::default()) {
            Ok(result) => {
                let total_solids: usize = result
                    .document
                    .products
                    .iter()
                    .flat_map(|p| &p.shapes)
                    .map(|s| s.solids.len())
                    .sum();
                println!("  STEP re-import OK: {total_solids} solid(s)");
            }
            Err(e) => eprintln!("  STEP re-import FAILED: {e}"),
        },
        Err(e) => eprintln!("  Could not read STEP file for validation: {e}"),
    }
}

/// Render a BRep from four camera angles, saving PNGs to `examples/output/`.
///
/// Produces `<base>_top.png`, `<base>_front.png`, `<base>_low.png`,
/// and `<base>_inside.png`.
#[allow(dead_code)]
pub fn render_multi_angle(brep: &BRep, base: &str) {
    use glam::Vec3;
    use popkern_viz::camera::Camera;

    let output_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("examples")
        .join("output");
    std::fs::create_dir_all(&output_dir).expect("create output dir");

    let base_camera = Camera {
        fov_y: std::f32::consts::FRAC_PI_4,
        aspect: 4.0 / 3.0,
        near: 0.01,
        far: 100.0,
        // eye/target/up are overwritten per-angle below
        eye: Vec3::ZERO,
        target: Vec3::ZERO,
        up: Vec3::Y,
    };

    let angles: &[(&str, Vec3, Vec3, Vec3)] = &[
        (
            "top",
            Vec3::new(0.5, 0.5, 4.0),
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(0.0, 1.0, 0.0),
        ),
        (
            "front",
            Vec3::new(0.5, -3.0, 0.5),
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(0.0, 0.0, 1.0),
        ),
        (
            "low",
            Vec3::new(2.0, -1.0, -0.5),
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(0.0, 0.0, 1.0),
        ),
        (
            "inside",
            Vec3::new(1.5, 1.5, 2.0),
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(0.0, 0.0, 1.0),
        ),
    ];

    for &(label, eye, target, up) in angles {
        let camera = Camera {
            eye,
            target,
            up,
            ..base_camera.clone()
        };
        let filename = format!("{base}_{label}.png");
        let path = output_dir.join(&filename);
        match popkern_viz::headless::render_to_png(brep, &camera, 800, 600, &path, true) {
            Ok(()) => println!("  Rendered: {}", path.display()),
            Err(e) => eprintln!("  Render {filename} failed: {e}"),
        }
    }
}

/// Try to render a BRep to PNG via the headless renderer.
/// Returns Ok if successful, or prints a warning and returns Err if GPU is unavailable.
pub fn try_render(brep: &BRep, filename: &str) {
    let output_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("examples")
        .join("output");
    std::fs::create_dir_all(&output_dir).expect("create output dir");
    let path = output_dir.join(filename);

    let camera = popkern_viz::camera::Camera::default_iso();
    match popkern_viz::headless::render_to_png(brep, &camera, 800, 600, &path, true) {
        Ok(()) => println!("  Rendered: {}", path.display()),
        Err(e) => eprintln!("  PNG render failed: {e}"),
    }
}
