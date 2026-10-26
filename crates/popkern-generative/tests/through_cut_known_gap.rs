//! PK-957 (known gap): exact through-cuts across a faceted design body.
//! The crisp block_two_bores design (`[generative]` fixture, 0.5 mm
//! cells) as a body, thresholded as the feature makes it: an exact boss
//! fused on its kept top and a hole drilled into the boss and the kept
//! plate serve (PK-956 P2, measured 6 of 6), but a hole cut through the
//! organic region declined at about half the positions in popkern-bool
//! (stranded half-claimed records, open loops, trims millimetres off).
//! The design takes minutes to solve, so the surface is cached in the
//! target's temporary directory.
//!
//! `cargo test --release -p popkern-generative --test through_cut_known_gap -- --ignored`

use std::path::PathBuf;

use popkern_analysis::Provenance;
use popkern_analysis::structural::from_setup::solve_generative_setup;
use popkern_analysis::structural::generative::body::design_model;
use popkern_analysis_brep::{FeatureTolerances, model_from_breps, resolve_faces};
use popkern_analysis_setup::parse_setup;
use popkern_core::tolerance::Tolerance;
use popkern_generative::facet::{FacetOptions, facet};
use popkern_io::{StepImportOptions, import_step};
use popkern_mesh::TessellationOptions;

fn repo(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

type Surface = (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<u32>);

fn model() -> popkern_analysis::AnalysisModel {
    let bytes = std::fs::read(repo("tests/fixtures/analysis/block_two_bores.step"))
        .unwrap_or_else(|e| panic!("{e}"));
    let imported =
        import_step(&bytes, &StepImportOptions::default()).unwrap_or_else(|e| panic!("{e}"));
    let breps: Vec<_> = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter().flat_map(|s| s.solids.iter().cloned()))
        .collect();
    model_from_breps(
        &breps,
        &TessellationOptions::display(),
        &FeatureTolerances::default(),
        Provenance::Step {
            file_sha256: String::new(),
        },
    )
    .unwrap_or_else(|e| panic!("{e}"))
}

/// The probe's faceting at the default tolerance.
fn plain_facet(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    owners: &[u32],
) -> popkern_generative::facet::FacetedBody {
    let tolerance = Tolerance::default();
    facet(
        vertices,
        triangles,
        owners,
        &FacetOptions::for_tolerance(&tolerance),
        tolerance,
    )
    .unwrap_or_else(|e| panic!("{e}"))
}

fn solve() -> Surface {
    let model = model();
    let text = std::fs::read_to_string(repo(
        "tests/fixtures/analysis/block_two_bores.generative.toml",
    ))
    .unwrap_or_else(|e| panic!("{e}"));
    let setup = parse_setup(&text).unwrap_or_else(|e| panic!("{e}"));
    let resolution =
        resolve_faces(&model, &setup.faces, &setup.tolerance).unwrap_or_else(|e| panic!("{e}"));
    let (design, _) =
        solve_generative_setup(&model, &setup, &resolution.faces).unwrap_or_else(|e| panic!("{e}"));
    let surface = popkern_generative::surface::body_surface(&design).mesh;
    let body = design_model(&model, &surface, design.cell_size);
    let owners = body.solids[0].triangle_face.iter().map(|f| f.0).collect();
    (surface.vertices, surface.triangles, owners)
}

fn cached() -> Surface {
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("block_two_bores.design.txt");
    if let Ok(text) = std::fs::read_to_string(&path) {
        let mut lines = text.lines();
        let mut next = || lines.next().unwrap_or_default();
        let nv: usize = next().parse().unwrap_or(0);
        let vertices = (0..nv)
            .map(|_| {
                let v: Vec<f64> = next()
                    .split(' ')
                    .map(|x| f64::from_bits(u64::from_str_radix(x, 16).unwrap_or(0)))
                    .collect();
                [v[0], v[1], v[2]]
            })
            .collect();
        let nt: usize = next().parse().unwrap_or(0);
        let mut triangles = Vec::new();
        let mut owners = Vec::new();
        for _ in 0..nt {
            let v: Vec<u32> = next().split(' ').map(|x| x.parse().unwrap_or(0)).collect();
            triangles.push([v[0], v[1], v[2]]);
            owners.push(v[3]);
        }
        return (vertices, triangles, owners);
    }
    let (vertices, triangles, owners) = solve();
    let mut text = format!("{}\n", vertices.len());
    for v in &vertices {
        text += &format!(
            "{:x} {:x} {:x}\n",
            v[0].to_bits(),
            v[1].to_bits(),
            v[2].to_bits()
        );
    }
    text += &format!("{}\n", triangles.len());
    for (t, o) in triangles.iter().zip(&owners) {
        text += &format!("{} {} {} {o}\n", t[0], t[1], t[2]);
    }
    std::fs::write(&path, text).unwrap_or_else(|e| panic!("{e}"));
    (vertices, triangles, owners)
}

fn faceted_part() -> (popkern_part::Part, [f64; 3], [f64; 3]) {
    let (vertices, triangles, owners) = cached();
    let tolerance = Tolerance::default();
    let body = plain_facet(&vertices, &triangles, &owners);
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in &vertices {
        for a in 0..3 {
            lo[a] = lo[a].min(v[a]);
            hi[a] = hi[a].max(v[a]);
        }
    }
    let part = popkern_part::Part::from_brep(body.brep, tolerance, "design")
        .unwrap_or_else(|e| panic!("{e}"));
    (part, lo, hi)
}

fn cylinder_at(x: f64, y: f64, z: f64, r: f64, h: f64) -> popkern_part::Part {
    use popkern_core::geometry::Transform3;
    let tolerance = Tolerance::default();
    let c = popkern_core::topology::BRep::cylinder(r, h, tolerance)
        .and_then(|b| b.transform(&Transform3::translation(x, y, z)))
        .unwrap_or_else(|e| panic!("{e}"));
    popkern_part::Part::from_brep(c, tolerance, "cylinder").unwrap_or_else(|e| panic!("{e}"))
}

#[test]
#[ignore = "known gap: PK-957 exact through-cuts across a faceted design body decline (solves the fixture's design first: minutes)"]
fn exact_through_cuts_across_the_design_serve() {
    let (part, lo, hi) = faceted_part();
    let mut declined = Vec::new();
    // Clear of the bores, through the kept top plate and the organic
    // region below it.
    for (dx, dy) in [
        (0.0, 2.6),
        (0.0, -2.6),
        (2.5, 3.3),
        (-2.7, -3.4),
        (0.13, 2.17),
    ] {
        let (cx, cy) = (0.5 * (lo[0] + hi[0]) + dx, 0.5 * (lo[1] + hi[1]) + dy);
        let drill = cylinder_at(cx, cy, lo[2] - 1.0, 1.2, hi[2] - lo[2] + 6.0);
        if let Err(e) = part.subtract(&drill) {
            declined.push(format!("({dx}, {dy}): {e}"));
        }
    }
    assert!(declined.is_empty(), "{}", declined.join("\n"));
}
