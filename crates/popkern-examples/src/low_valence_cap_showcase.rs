//! Low-valence corner caps — the n = 3 and n = 4 story.
//!
//! Valence tracks corner frequency in reverse. `n = 3` is the cube
//! corner and `n = 4` the pyramid apex or cross junction: between
//! them they account for most corners in real mechanical parts,
//! while `n ≥ 5` is comparatively rare. Yet the PnS bi-4 cap — the
//! kernel's best n-sided cap construction, G¹ across its sector
//! seams by construction — was dispatched only at `n ≥ 5` for most
//! of its life:
//!
//! - `n = 4` was believed structurally impossible, because σ's
//!   `σ_{20,x}` carries a `1/cos(2π/n)` factor that blows up at
//!   `n = 4`. It is not impossible: regrouped as
//!   `d₁ + (4/(3c))·(cos(α/2)·d₂ − d₁)`, the residue vanishes
//!   *bit-exactly* under the paper's Appendix B scalings, so the
//!   pole is removable and the jet is finite.
//! - `n = 3` was simply never wired into dispatch, though the
//!   adapter always supported it.
//!
//! This example measures what the low-valence corners actually ship,
//! writes the continuity numbers as a report, and renders the caps
//! so the result can be inspected rather than taken on trust.
//!
//! Outputs:
//!   examples/output/low_valence_cap_n3.png / .step
//!   examples/output/low_valence_cap_n4.png / .step
//!   examples/output/low_valence_cap_report.txt

use std::fmt::Write as _;
use std::path::Path;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_fillet::cap_seam_continuity::{measure_cap_seam_continuity, seam_normal_profile};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn main() {
    if let Err(e) = run() {
        eprintln!("low_valence_cap_showcase failed: {e}");
        std::process::exit(1);
    }
}

/// One corner's measured outcome.
struct CornerReport {
    valence: usize,
    radii: Vec<f64>,
    methods: String,
    achieved: String,
    cap_faces: usize,
    seams: usize,
    /// False when the cap has no internal seam at all — which must
    /// never be read as "perfectly continuous".
    measured: bool,
    worst_g0: f64,
    worst_normal_angle_deg: f64,
    /// Where along the seam the worst jump sits. Near 0.05 or 0.95
    /// points at the degenerate pole; mid-range means a real crease.
    worst_at_fraction: Option<f64>,
    mid_seam_angle_deg: f64,
}

fn run() -> Result<(), String> {
    let out_dir = Path::new("examples/output");
    std::fs::create_dir_all(out_dir).map_err(|e| format!("create output dir: {e}"))?;

    let mut report = String::new();
    writeln!(report, "Low-valence corner cap report").ok();
    writeln!(report, "=============================").ok();
    writeln!(
        report,
        "\nContinuity is measured across the cap's OWN internal seams\n\
         (cap sector vs adjacent cap sector), sampled away from the\n\
         degenerate pole. Acceptance in the exact-boundary route gates\n\
         cap-vs-SUPPORT coplanarity only, so this is the number that\n\
         route never checks.\n"
    )
    .ok();

    for (valence, radii) in [
        (3_usize, vec![0.05, 0.05, 0.05]),
        (3, vec![0.04, 0.06, 0.05]),
        (4, vec![0.05, 0.05, 0.05, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let (brep, corner) = build_apex(valence, &radii)?;
        writeln!(
            report,
            "n={} radii={:?}\n  route            : {}\n  achieved         : {}\n  cap faces        : {}\n  measured seams   : {} (measurable={})\n  worst G0 gap     : {:.3e}\n  worst normal jump: {:.1} deg{}\n  mid-seam jump    : {:.1} deg",
            corner.valence,
            corner.radii,
            corner.methods,
            corner.achieved,
            corner.cap_faces,
            corner.seams,
            corner.measured,
            corner.worst_g0,
            corner.worst_normal_angle_deg,
            corner
                .worst_at_fraction
                .map(|f| format!("  (at f={f:.2} along the seam)"))
                .unwrap_or_default(),
            corner.mid_seam_angle_deg,
        )
        .ok();

        let uniform = radii.windows(2).all(|w| (w[0] - w[1]).abs() < 1e-12);
        let tag = if uniform { "uniform" } else { "mixed" };
        let png = out_dir.join(format!("low_valence_cap_n{valence}_{tag}.png"));
        let r_max = radii.iter().copied().fold(0.0_f64, f64::max);
        render(&brep, &png, 3.0, r_max)?;
        writeln!(report, "  render           : {}", png.display()).ok();
        let step = out_dir.join(format!("low_valence_cap_n{valence}_{tag}.step"));
        match export_step(&brep, &step, valence) {
            Ok(()) => {
                writeln!(report, "  step             : {}", step.display()).ok();
            }
            Err(e) => {
                writeln!(report, "  step             : declined ({e})").ok();
            }
        }
        writeln!(report).ok();
    }

    let report_path = out_dir.join("low_valence_cap_report.txt");
    std::fs::write(&report_path, &report).map_err(|e| format!("write report: {e}"))?;
    print!("{report}");
    println!("report written to {}", report_path.display());
    Ok(())
}

/// Fillet the apex of an n-gonal pyramid and measure the cap.
fn build_apex(valence: usize, radii: &[f64]) -> Result<(BRep, CornerReport), String> {
    let tol = Tolerance::default();
    let pyramid =
        BRep::n_gonal_pyramid(valence, 1.0, 3.0, tol).map_err(|e| format!("pyramid: {e}"))?;
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
    if apex_edges.len() != valence {
        return Err(format!(
            "expected {valence} apex edges, found {}",
            apex_edges.len()
        ));
    }

    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .zip(radii)
        .map(|(&edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    let result =
        fillet_edges(&pyramid, &specs, &opts, tol, None).map_err(|e| format!("fillet: {e:?}"))?;

    let mut methods: Vec<String> = Vec::new();
    let mut achieved: Vec<String> = Vec::new();
    for patch in &result.corner_patches {
        let m = format!("{:?}", patch.method);
        if !methods.contains(&m) {
            methods.push(m);
        }
        let a = patch.achieved_continuity.as_str().to_string();
        if !achieved.contains(&a) {
            achieved.push(a);
        }
    }

    let cap_faces = cap_faces(&result.brep, &result.corner_patches);
    // Measure through the kernel's own cap-seam continuity API — the
    // same code the gates use — rather than re-deriving normals here
    // by point projection, which picks the wrong chart branch on a
    // folded patch and saturates at a meaningless 90 degrees.
    let seam = measure_cap_seam_continuity(&result.brep, &cap_faces, 100);
    let mid_angle = mid_seam_angle(&result.brep, &cap_faces);

    Ok((
        result.brep,
        CornerReport {
            valence,
            radii: radii.to_vec(),
            methods: if methods.is_empty() {
                "<no corner patches recorded>".to_string()
            } else {
                methods.join(", ")
            },
            achieved: if achieved.is_empty() {
                "-".to_string()
            } else {
                achieved.join(", ")
            },
            cap_faces: cap_faces.len(),
            seams: seam.seam_count,
            measured: seam.measured(),
            worst_g0: seam.worst_g0_gap,
            worst_normal_angle_deg: seam.worst_normal_angle_deg,
            worst_at_fraction: seam.worst_g1_fraction,
            mid_seam_angle_deg: mid_angle,
        },
    ))
}

fn cap_faces(brep: &BRep, patches: &[popkern_fillet::types::CornerPatchInfo]) -> Vec<FaceHandle> {
    let mut faces: Vec<FaceHandle> = Vec::new();
    for patch in patches {
        if brep.faces.contains_key(patch.face) && !faces.contains(&patch.face) {
            faces.push(patch.face);
        }
    }
    faces
}

/// Worst normal jump at the MIDDLE of each seam.
///
/// The distinction that matters: a jump concentrated at the seam ends
/// is a degenerate-pole artifact, while one present mid-seam is a
/// genuine crease in the shipped surface.
fn mid_seam_angle(brep: &BRep, faces: &[FaceHandle]) -> f64 {
    let mut worst = 0.0_f64;
    for (i, &fa) in faces.iter().enumerate() {
        for &fb in &faces[i + 1..] {
            if let Some(profile) = seam_normal_profile(brep, fa, fb, &[0.45, 0.50, 0.55]) {
                for (_f, angle) in profile {
                    worst = worst.max(angle);
                }
            }
        }
    }
    worst
}

fn export_step(brep: &BRep, path: &Path, valence: usize) -> Result<(), String> {
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: format!("low_valence_cap_n{valence}"),
            id: format!("LOW-VALENCE-CAP-N{valence}"),
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

/// Frame tightly on the apex cap.
///
/// The pyramid is 3 units tall and the fillet radius is ~0.05, so a
/// default view shows a speck. The camera is placed a few cap-radii
/// from the apex, looking slightly down, which is the only framing at
/// which a crease in the cap is actually visible.
fn render(brep: &BRep, path: &Path, apex_z: f64, radius: f64) -> Result<(), String> {
    let mut camera = popkern_viz::camera::Camera::default_iso();
    camera.aspect = 1.0;
    let r = radius as f32;
    let z = apex_z as f32;
    // Look straight down the pyramid axis, from just above the cap.
    // The cap spans roughly 2r, so a 3r standoff with a 45 degree
    // field fills the frame with it.
    // The filleted apex sits about 1.5r below the original point.
    camera.target = glam::Vec3::new(0.0, 0.0, z - r * 1.6);
    camera.eye = glam::Vec3::new(r * 1.7, -r * 1.7, z + r * 1.2);
    camera.up = glam::Vec3::Z;
    camera.fov_y = 45.0_f32.to_radians();
    camera.near = r * 0.05;
    popkern_viz::headless::render_to_png(brep, &camera, 900, 900, path, true)
        .map_err(|e| format!("render {}: {e}", path.display()))
}

/// Surface-kind census, printed so the reader can see which
/// construction actually shipped rather than inferring it.
#[allow(dead_code)]
fn surface_kinds(brep: &BRep, faces: &[FaceHandle]) -> Vec<&'static str> {
    faces
        .iter()
        .filter_map(|&f| brep.faces.get(f)?.surface.as_ref())
        .map(|s| match &s.geometry {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cylinder",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
            SurfaceGeometry::Gregory(_) => "Gregory",
        })
        .collect()
}
