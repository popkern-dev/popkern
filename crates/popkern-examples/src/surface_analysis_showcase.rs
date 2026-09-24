//! Surface-analysis showcase (PK-520): the pentagon-apex G2 flagship
//! rendered in the three inspection views — shaded, reflection-line
//! zebra, and mean-curvature false-color — from a camera looking down
//! at the setback vertex cap. The zebra stripes flowing unbroken from
//! the fillet strips across the rim seams into the cap are the visual
//! counterpart of the PK-519 rim G2 gate.
//!
//! Outputs:
//!   examples/output/surface_analysis_shaded.png
//!   examples/output/surface_analysis_zebra.png
//!   examples/output/surface_analysis_gaussian.png
//!   examples/output/surface_analysis_mean.png

use std::path::Path;

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::CurvatureKind;
use popkern_viz::camera::Camera;
use popkern_viz::surface_analysis::{SurfaceAnalysisMode, render_surface_analysis_to_png};

fn main() {
    println!("=== Surface Analysis Showcase ===");
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
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
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (origin == apex || twin == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.16),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let result = fillet_edges(&pyramid, &specs, &options, tol, None).expect("fillet");
    let brep = result.brep;
    println!(
        "  Flagship: {} faces, {} corner patches",
        brep.faces.len(),
        result.corner_patches.len()
    );

    // Per-face curvature census: surface kind + mean-curvature range
    // over each face's mesh vertices (the numbers behind the map).
    {
        let mesh = popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::default())
            .expect("mesh");
        let field =
            popkern_measure::curvature_field(&brep, &mesh, CurvatureKind::Mean).expect("field");
        let mut cursor = 0usize;
        for fm in &mesh.face_meshes {
            let n = fm.vertices.len();
            let slice = &field.per_vertex[cursor..cursor + n];
            cursor += n;
            let kind = brep
                .faces
                .get(fm.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.kind_name())
                .unwrap_or("none");
            let (lo, hi) = slice
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
                    (a.min(v), b.max(v))
                });
            // Direct per-face evaluation (own surface, own uv) as the
            // reference the combined-mesh field must match.
            let (dlo, dhi) = brep
                .faces
                .get(fm.face)
                .and_then(|f| f.surface.as_ref())
                .map(|sof| {
                    fm.vertices
                        .iter()
                        .filter_map(|v| {
                            popkern_measure::curvature_at(
                                &sof.geometry,
                                v.uv.0,
                                v.uv.1,
                                CurvatureKind::Mean,
                            )
                        })
                        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| {
                            (a.min(v), b.max(v))
                        })
                })
                .unwrap_or((f64::NAN, f64::NAN));
            println!(
                "  face {:?} {kind}: H in [{lo:.3}, {hi:.3}] (direct [{dlo:.3}, {dhi:.3}]) over {n} vertices",
                fm.face
            );
        }
    }

    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/output");
    std::fs::create_dir_all(&out_dir).expect("output dir");

    // Close-up on the setback vertex cap: the rim seams and their
    // stripe/curvature continuity are the subject.
    let camera = Camera {
        eye: glam::Vec3::new(0.75, -0.62, 1.75),
        target: glam::Vec3::new(0.0, 0.0, 1.02),
        up: glam::Vec3::Z,
        fov_y: std::f32::consts::FRAC_PI_4,
        aspect: 4.0 / 3.0,
        near: 0.01,
        far: 100.0,
    };

    let shaded = out_dir.join("surface_analysis_shaded.png");
    popkern_viz::headless::render_to_png_software(&brep, &camera, 800, 600, &shaded, false)
        .expect("shaded render");
    println!("  wrote {}", shaded.display());

    for (name, mode) in [
        (
            "surface_analysis_zebra.png",
            SurfaceAnalysisMode::Zebra { stripe_count: 14 },
        ),
        (
            "surface_analysis_gaussian.png",
            SurfaceAnalysisMode::Curvature {
                kind: CurvatureKind::Gaussian,
            },
        ),
        (
            "surface_analysis_mean.png",
            SurfaceAnalysisMode::Curvature {
                kind: CurvatureKind::Mean,
            },
        ),
    ] {
        let path = out_dir.join(name);
        render_surface_analysis_to_png(&brep, &camera, 800, 600, mode, &path)
            .expect("analysis render");
        println!("  wrote {}", path.display());
    }
}
