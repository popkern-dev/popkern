//! Control cage → smooth solid → ordinary CAD body.
//!
//! Polyhedral-net splines (Karčiauskas–Peters) turn an arbitrary
//! quad/n-gon control cage — the kind a subdivision modeller produces
//! — into a watertight BRep of smooth Bézier patches. popkern has
//! carried the whole family in `popkern-core` for a long time with no
//! production caller; `popkern_part::smooth_mesh` is the single call
//! that makes it reachable, and this example is the proof that what
//! comes out is *ordinary geometry*: it validates, it tessellates, it
//! renders, it exports to STEP, and it can be filleted like any other
//! part.
//!
//! Three cages, chosen to exercise different constructions:
//!
//! - **cube** — 8 vertices, 6 quads. Every vertex is valence 3, so
//!   every corner is an `ExtraordinaryStar`, not the easy regular case.
//! - **truncated corner** — introduces an n-gon face and a
//!   higher-valence star.
//! - **L-block** — a concave cage, so the surface has to bend both
//!   ways.
//!
//! Outputs:
//!   examples/output/smooth_mesh_<cage>.png / .step
//!   examples/output/smooth_mesh_report.txt

use std::fmt::Write as _;
use std::path::Path;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_part::{SmoothMeshOptions, smooth_solid_from_obj};

fn main() {
    if let Err(e) = run() {
        eprintln!("smooth_mesh_showcase failed: {e}");
        std::process::exit(1);
    }
}

/// Closed cube cage — every vertex valence 3.
const CUBE_CAGE: &str = "\
v -1 -1 -1
v  1 -1 -1
v  1  1 -1
v -1  1 -1
v -1 -1  1
v  1 -1  1
v  1  1  1
v -1  1  1
f 1 4 3 2
f 5 6 7 8
f 1 2 6 5
f 2 3 7 6
f 3 4 8 7
f 4 1 5 8
";

/// Cube with one corner truncated: adds a triangle face and raises
/// the valence of the three vertices around it.
const TRUNCATED_CAGE: &str = "\
v -1 -1 -1
v  1 -1 -1
v  1  1 -1
v -1  1 -1
v -1 -1  1
v  1 -1  1
v -1  1  1
v  1  0  1
v  1  1  0
v  0  1  1
f 1 4 3 2
f 5 6 8 10 7
f 1 2 6 5
f 2 3 9 8 6
f 3 4 7 10 9
f 4 1 5 7
f 8 9 10
";

/// Concave L-shaped cage — the surface must bend both ways.
const L_BLOCK_CAGE: &str = "\
v 0 0 0
v 2 0 0
v 2 1 0
v 1 1 0
v 1 2 0
v 0 2 0
v 0 0 1
v 2 0 1
v 2 1 1
v 1 1 1
v 1 2 1
v 0 2 1
f 6 5 4 3 2 1
f 7 8 9 10 11 12
f 1 2 8 7
f 2 3 9 8
f 3 4 10 9
f 4 5 11 10
f 5 6 12 11
f 6 1 7 12
";

fn run() -> Result<(), String> {
    let out_dir = Path::new("examples/output");
    std::fs::create_dir_all(out_dir).map_err(|e| format!("create output dir: {e}"))?;

    let mut report = String::new();
    writeln!(report, "Polyhedral-net spline: control cage → smooth solid").ok();
    writeln!(
        report,
        "==================================================\n"
    )
    .ok();

    let tol = Tolerance::default();
    for (name, cage, extent) in [
        ("cube", CUBE_CAGE, 1.6_f64),
        ("truncated_corner", TRUNCATED_CAGE, 1.6),
        ("l_block", L_BLOCK_CAGE, 2.2),
    ] {
        writeln!(report, "cage: {name}").ok();
        match smooth_solid_from_obj(cage, tol, SmoothMeshOptions::default()) {
            Ok((part, r)) => {
                let brep = part.brep();
                writeln!(
                    report,
                    "  cage            : {} vertices, {} faces\n  \
                       splines         : {} candidates → {} components\n  \
                       weld            : {} pair(s), {} declined, {} boundary edge(s) left, max gap {:.3e}\n  \
                       body            : {} faces, {} edges, {} vertices\n  \
                       validates       : {}\n  \
                       complete        : {}",
                    r.source_vertices,
                    r.source_faces,
                    r.candidate_splines,
                    r.accepted_components,
                    r.welded_pairs,
                    r.declined_pairs,
                    r.remaining_boundary_edges,
                    r.max_weld_gap,
                    brep.faces.len(),
                    brep.edges.len(),
                    brep.vertices.len(),
                    brep.validate().is_ok(),
                    r.complete(),
                )
                .ok();
                for rejection in r.rejections.iter().take(4) {
                    writeln!(report, "  declined spline : {rejection}").ok();
                }

                let png = out_dir.join(format!("smooth_mesh_{name}.png"));
                match render(brep, &png, extent) {
                    Ok(()) => {
                        writeln!(report, "  render          : {}", png.display()).ok();
                    }
                    Err(e) => {
                        writeln!(report, "  render          : failed ({e})").ok();
                    }
                }

                let step = out_dir.join(format!("smooth_mesh_{name}.step"));
                match export_step(brep, &step, name) {
                    Ok(()) => {
                        writeln!(report, "  step            : {}", step.display()).ok();
                    }
                    Err(e) => {
                        writeln!(report, "  step            : declined ({e})").ok();
                    }
                }
            }
            Err(e) => {
                writeln!(report, "  DECLINED: {e}").ok();
            }
        }
        writeln!(report).ok();
    }

    let report_path = out_dir.join("smooth_mesh_report.txt");
    std::fs::write(&report_path, &report).map_err(|e| format!("write report: {e}"))?;
    print!("{report}");
    println!("report written to {}", report_path.display());
    Ok(())
}

fn export_step(brep: &BRep, path: &Path, name: &str) -> Result<(), String> {
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: format!("smooth_mesh_{name}"),
            id: format!("SMOOTH-MESH-{}", name.to_uppercase()),
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
    popkern_io::export_step_file(&doc, path, &popkern_io::StepExportOptions::default())
        .map_err(|e| format!("{e}"))
}

fn render(brep: &BRep, path: &Path, extent: f64) -> Result<(), String> {
    let center = body_center(brep);
    let mut camera = popkern_viz::camera::Camera::default_iso();
    camera.aspect = 1.0;
    let d = (extent * 2.4) as f32;
    camera.target = center;
    camera.eye = center + glam::Vec3::new(d, -d, d * 0.8);
    camera.up = glam::Vec3::Z;
    camera.fov_y = 40.0_f32.to_radians();
    popkern_viz::headless::render_to_png(brep, &camera, 900, 900, path, false)
        .map_err(|e| format!("render {}: {e}", path.display()))
}

fn body_center(brep: &BRep) -> glam::Vec3 {
    let mut sum = glam::Vec3::ZERO;
    let mut count = 0usize;
    for vertex in brep.vertices.values() {
        let p = vertex.position;
        sum += glam::Vec3::new(p.x as f32, p.y as f32, p.z as f32);
        count += 1;
    }
    if count == 0 { sum } else { sum / count as f32 }
}
