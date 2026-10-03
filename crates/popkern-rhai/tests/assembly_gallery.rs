//! PK-629 acceptance: the gallery's multi-part models rewritten as
//! assemblies build the SAME geometry as the frozen single-tree
//! scripts they replace (`tests/fixtures/single_model_*.rhai`, the gallery
//! before the rewrite), and build each repeated part once.
//!
//! Identity, per single-model body: exactly one flattened occurrence with
//! the same face, edge and vertex counts and the same surface kinds;
//! every vertex that is not on a seam edge coincident (to 1e-9 × model
//! size) with one of the other body, both ways; and the same volume to
//! the measure's tessellation accuracy. Seam vertices are excluded
//! because they are a parametrisation, not the geometry: a bore circle
//! cut at its centre and one turned into place seam at different
//! angles. Where the rewrite keeps the single-model operation order (the
//! bearing's balls, the fastener kit) every vertex position is
//! bit-identical. A planet turned one degree off its phase fails
//! (the negative control, measured when this test was written).

use std::collections::BTreeMap;

use popkern_assembly::{Assembly, Document};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_rhai::ScriptHost;
use popkern_script::Model;

fn root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

fn single_model(stem: &str) -> Model {
    let source =
        std::fs::read_to_string(root().join(format!("tests/fixtures/single_model_{stem}.rhai")))
            .expect("single-model fixture");
    let mut outcome = ScriptHost::new()
        .run(&source, stem)
        .expect("single-model script runs");
    match outcome.take_document() {
        Some(Document::Part(model)) => *model,
        _ => panic!("{stem}: the single-model script authors one model"),
    }
}

fn rewritten(stem: &str) -> Assembly {
    let source = std::fs::read_to_string(root().join(format!("examples/{stem}.rhai")))
        .expect("gallery script");
    let mut outcome = ScriptHost::new().run(&source, stem).expect("assembly runs");
    match outcome.take_document() {
        Some(Document::Assembly(assembly)) => *assembly,
        _ => panic!("{stem}: the gallery script authors an assembly"),
    }
}

struct Signature {
    counts: (usize, usize, usize),
    kinds: BTreeMap<String, usize>,
    vertices: Vec<[f64; 3]>,
    centroid: [f64; 3],
}

fn signature(brep: &BRep) -> Signature {
    let mut kinds = BTreeMap::new();
    for face in brep.faces.values() {
        let kind = face
            .surface
            .as_ref()
            .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)))
            .unwrap_or_default();
        *kinds.entry(kind).or_insert(0) += 1;
    }
    let vertices: Vec<[f64; 3]> = brep
        .vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect();
    let n = vertices.len().max(1) as f64;
    let mut centroid = [0.0; 3];
    for v in &vertices {
        for k in 0..3 {
            centroid[k] += v[k] / n;
        }
    }
    Signature {
        counts: (brep.faces.len(), brep.edges.len(), brep.vertices.len()),
        kinds,
        vertices,
        centroid,
    }
}

fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// The largest distance from a vertex of `a` to its nearest vertex of `b`.
fn covered(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter()
        .map(|p| {
            b.iter()
                .map(|q| distance(p, q))
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max)
}

/// The vertices of `brep` that are not on a seam edge (an edge whose
/// two sides bound the same face).
fn feature_vertices(brep: &BRep) -> Vec<[f64; 3]> {
    let mut seam_vertices = std::collections::HashSet::new();
    for edge in brep.edges.values() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let face = |h: &popkern_core::topology::HalfEdge| brep.loops.get(h.loop_).map(|l| l.face);
        if face(he).is_some() && face(he) == face(twin) {
            seam_vertices.insert(he.origin);
            seam_vertices.insert(twin.origin);
        }
    }
    brep.vertices
        .iter()
        .filter(|(h, _)| !seam_vertices.contains(h))
        .map(|(_, v)| [v.position.x, v.position.y, v.position.z])
        .collect()
}

struct Comparison {
    bodies: usize,
    worst_vertex: f64,
    worst_volume: f64,
    bit_identical: bool,
}

fn compare(stem: &str) -> (Comparison, Model, Assembly) {
    let mut old = single_model(stem);
    let mut new = rewritten(stem);
    let old_bodies: Vec<BRep> = old
        .to_parts_all()
        .expect("single-model parts")
        .into_iter()
        .map(|p| p.brep().clone())
        .collect();
    let flat = new.flatten().expect("flatten");
    assert_eq!(
        flat.len(),
        old_bodies.len(),
        "[{stem}] one occurrence per single-model body"
    );
    let size = old_bodies
        .iter()
        .flat_map(|b| {
            b.vertices.values().map(|v| {
                v.position
                    .x
                    .abs()
                    .max(v.position.y.abs())
                    .max(v.position.z.abs())
            })
        })
        .fold(0.0f64, f64::max);
    let tol = 1e-9 * size;
    let mut used = vec![false; flat.len()];
    let mut worst_vertex: f64 = 0.0;
    let mut worst_volume: f64 = 0.0;
    let mut bit_identical = true;
    for old_body in &old_bodies {
        let want = signature(old_body);
        let found = flat
            .iter()
            .enumerate()
            .filter(|(i, f)| !used[*i] && signature(&f.brep).counts == want.counts)
            .min_by(|(_, x), (_, y)| {
                let dx = distance(&signature(&x.brep).centroid, &want.centroid);
                let dy = distance(&signature(&y.brep).centroid, &want.centroid);
                dx.total_cmp(&dy)
            });
        let Some((i, body)) = found else {
            panic!(
                "[{stem}] no occurrence matches a single-model body with counts {:?}",
                want.counts
            )
        };
        used[i] = true;
        let got = signature(&body.brep);
        assert_eq!(
            got.kinds, want.kinds,
            "[{stem}] {}: surface kinds",
            body.path
        );
        let there = covered(&want.vertices, &got.vertices);
        let back = covered(&got.vertices, &want.vertices);
        let (there, back) = if there <= tol && back <= tol {
            (there, back)
        } else {
            let (old_features, new_features) =
                (feature_vertices(old_body), feature_vertices(&body.brep));
            assert_eq!(
                old_features.len(),
                new_features.len(),
                "[{stem}] {}: feature vertices",
                body.path
            );
            (
                covered(&old_features, &new_features),
                covered(&new_features, &old_features),
            )
        };
        worst_vertex = worst_vertex.max(there).max(back);
        let mut a = want.vertices.clone();
        let mut b = got.vertices.clone();
        a.sort_by(|p, q| p.partial_cmp(q).expect("finite"));
        b.sort_by(|p, q| p.partial_cmp(q).expect("finite"));
        let same_bits = a == b;
        bit_identical &= same_bits;
        if !same_bits {
            // Equal bits are the same solid; otherwise measure it.
            let tolerance = Tolerance::default();
            let va = popkern_measure::body_volume(old_body, &tolerance).expect("volume");
            let vb = popkern_measure::body_volume(&body.brep, &tolerance).expect("volume");
            worst_volume = worst_volume.max(((va - vb) / va).abs());
        }
    }
    assert!(
        worst_vertex <= tol,
        "[{stem}] worst vertex deviation {worst_vertex:e} > {tol:e}"
    );
    (
        Comparison {
            bodies: old_bodies.len(),
            worst_vertex,
            worst_volume,
            bit_identical,
        },
        old,
        new,
    )
}

/// Each part definition's occurrence count, by name.
fn occurrences(asm: &Assembly) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for row in asm.leaves() {
        let name = asm
            .definition(row.definition)
            .expect("definition")
            .name
            .clone();
        *counts.entry(name).or_insert(0) += 1;
    }
    counts
}

fn features(asm: &mut Assembly) -> usize {
    asm.part_models_mut()
        .into_iter()
        .map(|(_, m)| m.features().len())
        .sum()
}

#[test]
fn the_bearing_assembly_is_the_single_model_bearing_bit_for_bit() {
    let (c, old, mut new) = compare("ball_bearing");
    assert_eq!(c.bodies, 10);
    assert!(c.bit_identical, "same operation order, same bits");
    assert!(c.worst_volume < 1e-12, "{:e}", c.worst_volume);
    assert_eq!(occurrences(&new)["ball"], 8, "one ball, eight occurrences");
    // Build count: eight balls (sketch, revolve, translate: 24
    // features) became one (sketch, revolve); each race is three.
    assert_eq!(old.features().len(), 3 + 3 + 24);
    assert_eq!(features(&mut new), 3 + 3 + 2);
}

#[test]
fn the_gearbox_assembly_is_the_single_model_gearbox() {
    let (c, old, mut new) = compare("planetary_gearbox");
    assert_eq!(c.bodies, 6);
    eprintln!(
        "gearbox: worst vertex {:e}, worst volume {:e} (relative), bit-identical {}",
        c.worst_vertex, c.worst_volume, c.bit_identical
    );
    assert!(c.worst_volume < 1e-6, "{:e}", c.worst_volume);
    assert_eq!(
        occurrences(&new)["planet"],
        3,
        "one planet, three occurrences"
    );
    // Build count: three planet outlines, extrudes and bore booleans
    // became one.
    let (single, now) = (old.features().len(), features(&mut new));
    assert_eq!(
        single - now,
        2 * 5,
        "two planets' outline sketch, extrude, bore sketch, bore extrude and bore boolean saved"
    );
}

#[test]
fn the_fastener_kit_assembly_is_the_single_model_kit_from_one_family() {
    let (c, _, new) = compare("fastener_kit");
    assert_eq!(c.bodies, 4);
    eprintln!(
        "fastener kit: worst vertex {:e}, worst volume {:e}",
        c.worst_vertex, c.worst_volume
    );
    assert!(c.worst_volume < 1e-6, "{:e}", c.worst_volume);
    let families: Vec<String> = new
        .definitions()
        .iter()
        .filter_map(|d| d.family.as_ref().map(|f| f.family.clone()))
        .collect();
    assert_eq!(
        families,
        ["screw", "screw", "screw"],
        "three screws, one family"
    );
}
