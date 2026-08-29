//! KP20 pyramid-apex showcase — the honest production pipeline.
//!
//! Builds a pentagonal pyramid, fillets its five apex edges with G²
//! sweeps through the production `fillet_edges` pipeline (whose
//! vertex blend is the KP20 §4.4 setback cap: machine-zero paper
//! gates, real trimmed topology, `validate()`d), tessellates the
//! blended solid with popkern-mesh, and exports OBJ + STEP.
//!
//! The exploratory probes that drove the campaign (fillet side
//! checks, tangency, curvature combs, chart verification) are
//! promoted to unit/integration tests in popkern-fillet; this
//! example is the clean production story.
//!
//! Outputs (picked up by tools/render-apex-showcase.sh):
//!   examples/output/kp20_pyramid_apex_whole.obj
//!   examples/output/kp20_pyramid_apex_blend_detail.obj
//!   examples/output/kp20_pyramid_apex_whole.step

use std::path::Path;

use popkern_core::geometry::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_mesh::{TessellationOptions, TriMesh, tessellate};

fn main() {
    if let Err(e) = run() {
        eprintln!("kp20_pyramid_apex_showcase failed: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).map_err(|e| format!("pyramid: {e}"))?;

    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| {
            a.1.position
                .z
                .partial_cmp(&b.1.position.z)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(h, _)| h)
        .ok_or("no vertices")?;
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (origin == apex || twin_origin == apex).then_some(eh)
        })
        .collect();
    println!("apex has {} incident edges", apex_edges.len());

    // Default radius r = 0.16 — the original flagship look. The
    // large-radius fan class declined on the fan-seam G1
    // postcondition until PK-094 landed the strip-movement-scaled
    // rim gate + cross-seam jet pinning (measured: g1 1.478e-1 →
    // 1.0e-3 on this exact fixture; the PK-123 serving boundary at
    // r = 0.06 is retired). Override with `--radius <r>` to probe.
    let radius: f64 = cli_value("--radius")
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.16);
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let continuity = match cli_value("--continuity").as_deref() {
        Some("G1") | Some("g1") => ContinuityLevel::G1,
        _ => ContinuityLevel::G2,
    };
    let options = FilletOptions {
        continuity,
        ..Default::default()
    };

    // The production pipeline: G² edge fillets + the KP20 setback
    // vertex blend (on by default for eligible apexes).
    let result = fillet_edges(&pyramid, &specs, &options, tol, None)
        .map_err(|e| format!("fillet_edges: {e:?}"))?;
    let brep = &result.brep;
    brep.validate()
        .map_err(|e| format!("blended BRep failed validation: {e}"))?;

    println!(
        "blended BRep: faces={} edges={} vertices={} (validate: ok)",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
    let gregory_faces = brep
        .faces
        .iter()
        .filter(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Gregory(_))
            )
        })
        .count();
    println!(
        "corner patches: {} (KP20 Gregory cap faces: {gregory_faces})",
        result.corner_patches.len()
    );
    // The showcase exists to exhibit the KP20 apex cap: a silent decline
    // to a mitred corner must fail loudly instead of exporting STEP with
    // no cap (PK-123).
    if result.corner_patches.is_empty() || gregory_faces == 0 {
        return Err(format!(
            "apex blend declined to a mitred corner (corner_patches={}, gregory_faces={gregory_faces});              the requested radius {radius} is outside the class the setback surgery serves",
            result.corner_patches.len()
        ));
    }
    if let Some(cp) = result.corner_patches.first() {
        println!(
            "  apex blend: n_sides={} method={:?} continuity={:?}",
            cp.n_sides, cp.method, cp.achieved_continuity
        );
    }
    for (_, w) in result.edge_warnings.iter().take(5) {
        println!("  warn: {w:?}");
    }

    // Tessellate through the kernel and export.
    let mesh = tessellate(brep, &TessellationOptions::default())
        .map_err(|e| format!("tessellation: {e}"))?;
    println!(
        "tessellated: {} vertices, {} triangles",
        mesh.vertices.len(),
        mesh.triangles.len()
    );

    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/output");
    std::fs::create_dir_all(&out_dir).map_err(|e| e.to_string())?;

    let whole_path = out_dir.join("kp20_pyramid_apex_whole.obj");
    write_grouped_obj(brep, &mesh, &whole_path, /*detail_only=*/ false)?;
    println!("wrote {}", whole_path.display());

    let detail_path = out_dir.join("kp20_pyramid_apex_blend_detail.obj");
    write_grouped_obj(brep, &mesh, &detail_path, /*detail_only=*/ true)?;
    println!("wrote {}", detail_path.display());

    let step_path = out_dir.join("kp20_pyramid_apex_whole.step");
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: "kp20_pyramid_apex".into(),
            id: "KP20-PYRAMID-APEX-001".into(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: popkern_io::RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata::default(),
        }],
        roots: vec![0],
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    match popkern_io::export_step_file(&doc, &step_path, &popkern_io::StepExportOptions::default())
    {
        Ok(()) => println!("wrote {}", step_path.display()),
        Err(e) => println!("STEP export declined: {e} (recorded: PK-094 tail)"),
    }
    Ok(())
}

/// Write the tessellation as an OBJ with one group per originating
/// face variant (fillet_N / paper_cap_N / face / other), so the
/// render script's blend-detail filter keeps working. `detail_only`
/// keeps just the fillet and cap groups.
fn write_grouped_obj(
    brep: &BRep,
    mesh: &TriMesh,
    path: &Path,
    detail_only: bool,
) -> Result<(), String> {
    use std::io::Write;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path).map_err(|e| e.to_string())?);
    writeln!(
        file,
        "# popkern kernel tessellation (fillet_edges + KP20 setback vertex blend)"
    )
    .map_err(|e| e.to_string())?;
    for v in &mesh.vertices {
        writeln!(file, "v {} {} {}", v.position.x, v.position.y, v.position.z)
            .map_err(|e| e.to_string())?;
        writeln!(file, "vn {} {} {}", v.normal.x, v.normal.y, v.normal.z)
            .map_err(|e| e.to_string())?;
    }
    let mut fillet_counter = 0usize;
    let mut cap_counter = 0usize;
    let mut group_of_face = std::collections::HashMap::new();
    let mut current_group = String::new();
    for (tri, &face) in mesh.triangles.iter().zip(&mesh.face_map) {
        let group = group_of_face
            .entry(face)
            .or_insert_with(|| {
                let kind = brep
                    .faces
                    .get(face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| &s.geometry);
                match kind {
                    Some(SurfaceGeometry::Gregory(_)) => {
                        cap_counter += 1;
                        format!("paper_cap_{}", cap_counter - 1)
                    }
                    Some(SurfaceGeometry::Nurbs(_)) | Some(SurfaceGeometry::Cylinder { .. }) => {
                        fillet_counter += 1;
                        format!("fillet_{}", fillet_counter - 1)
                    }
                    Some(SurfaceGeometry::Plane { .. }) => "face".to_string(),
                    _ => "other".to_string(),
                }
            })
            .clone();
        if detail_only && !(group.starts_with("fillet_") || group.starts_with("paper_cap_")) {
            continue;
        }
        if group != current_group {
            writeln!(file, "g {group}").map_err(|e| e.to_string())?;
            current_group = group;
        }
        writeln!(
            file,
            "f {}//{} {}//{} {}//{}",
            tri[0] + 1,
            tri[0] + 1,
            tri[1] + 1,
            tri[1] + 1,
            tri[2] + 1,
            tri[2] + 1
        )
        .map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// `--flag value` lookup over the process args.
fn cli_value(flag: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == flag)
        .and_then(|i| args.get(i + 1))
        .cloned()
}
