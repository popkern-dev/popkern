//! NOT COMPILED: copy to crates/popkern-blend/tests/zz_probe_corner_census.rs to run
//! (it is #[ignore]d; the ignore ratchet forbids committing it as a test).
//! Run: POPKERN_TRACE=corner-carve cargo test -p popkern-blend --test
//! zz_probe_corner_census -- --ignored --nocapture --exact census_pyramid_mixed
//! DIAGNOSTIC PROBE (re-run for PK-591..598; transient): which corners fall past
//! the KP20 setback cap today, and what serves them instead. Removed
//! once its findings are promoted into the PK-574 design doc.
//!
//! Run: POPKERN_TRACE=corner-carve cargo test -p popkern-blend --test
//! zz_probe_corner_census -- --ignored --nocapture --test-threads=1

use std::collections::{BTreeMap, BTreeSet};

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::Point3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

fn edges_at(body: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    edges
}

fn apex(body: &BRep) -> VertexHandle {
    body.vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex")
}

fn vertex_at(body: &BRep, at: Point3) -> VertexHandle {
    body.vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&at) < 1e-6)
        .map(|(h, _)| h)
        .expect("vertex")
}

fn run(name: &str, body: &BRep, specs: Vec<FilletSpec>, continuity: ContinuityLevel) {
    let tol = Tolerance::default();
    let options = FilletOptions {
        continuity,
        ..Default::default()
    };
    let pre: BTreeSet<FaceHandle> = body.faces.keys().collect();
    let _ = popkern_blend::construct_test_hooks::take_setback_topology_checkpoint();
    let t0 = std::time::Instant::now();
    let result = fillet_edges(body, &specs, &options, tol, None);
    let secs = t0.elapsed().as_secs_f64();
    let last = popkern_blend::construct_test_hooks::take_setback_topology_checkpoint();
    match result {
        Err(e) => eprintln!("CENSUS {name} {continuity:?}: FILLET ERROR {e:?} ({secs:.1}s)"),
        Ok(result) => {
            let mut per_vertex: BTreeMap<String, Vec<String>> = BTreeMap::new();
            for p in &result.corner_patches {
                per_vertex
                    .entry(format!("{:?}", p.vertex))
                    .or_default()
                    .push(format!(
                        "{:?}/{:?}/sb{}/n{}",
                        p.method, p.achieved_continuity, p.setback_cap as u8, p.n_sides
                    ));
            }
            let blend: Vec<FaceHandle> = result
                .brep
                .faces
                .keys()
                .filter(|f| !pre.contains(f))
                .collect();
            let ev = measure_blend_continuity_evidence(
                &result.brep,
                &blend,
                continuity,
                &ContinuityGates::commercial(),
                33,
            );
            eprintln!(
                "CENSUS {name} {continuity:?}: achieved {:?} g0 {:.2e} g1 {:.2e} g2 {:.2e} valid {} ({secs:.1}s)",
                ev.achieved,
                ev.worst_g0_gap,
                ev.worst_g1_normal,
                ev.worst_g2_residual,
                result.brep.validate().is_ok(),
            );
            for (v, list) in per_vertex {
                let mut counts: BTreeMap<String, usize> = BTreeMap::new();
                for s in list {
                    *counts.entry(s).or_default() += 1;
                }
                eprintln!("CENSUS   {v}: {counts:?}");
            }
            for d in &result.repair_report.apex_diagnostics {
                eprintln!("CENSUS   diag: {d:?}");
            }
            for (e, w) in &result.edge_warnings {
                eprintln!("CENSUS   warn {e:?}: {w:?}");
            }
            if let Some(last) = last {
                eprintln!("CENSUS   setback-last: {last}");
            }
        }
    }
}

fn pyramid(name: &str, n: usize, height: f64, radii: &[f64], continuity: ContinuityLevel) {
    let tol = Tolerance::default();
    let body = BRep::n_gonal_pyramid(n, 1.0, height, tol).expect("pyramid");
    let a = apex(&body);
    let specs = edges_at(&body, a)
        .into_iter()
        .zip(radii.iter().cycle())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    run(name, &body, specs, continuity);
}

#[test]
#[ignore = "diagnostic probe"]
fn census_pyramid_apexes() {
    for &n in &[3usize, 4, 5, 6, 7, 8] {
        for &h in &[0.25, 0.8, 1.2, 3.0] {
            for &r in &[0.05, 0.16, 0.3] {
                pyramid(
                    &format!("pyr{n}_h{h}_r{r}"),
                    n,
                    h,
                    &[r],
                    ContinuityLevel::G2,
                );
            }
        }
    }
}

#[test]
#[ignore = "diagnostic probe"]
fn census_pyramid_mixed() {
    for &n in &[3usize, 4, 5, 6, 8] {
        for &h in &[0.4, 0.8, 1.2] {
            let mixed = [0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06];
            let variable = [0.02, 0.04, 0.06, 0.04, 0.02, 0.05, 0.03, 0.06];
            let wide = [0.05, 0.15, 0.08, 0.2, 0.06, 0.12, 0.1, 0.18];
            pyramid(
                &format!("pyr{n}_h{h}_mixed"),
                n,
                h,
                &mixed,
                ContinuityLevel::G2,
            );
            pyramid(
                &format!("pyr{n}_h{h}_variable"),
                n,
                h,
                &variable,
                ContinuityLevel::G2,
            );
            pyramid(
                &format!("pyr{n}_h{h}_wide"),
                n,
                h,
                &wide,
                ContinuityLevel::G2,
            );
        }
    }
}

#[test]
#[ignore = "diagnostic probe"]
fn census_cube_and_saddle() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let c = vertex_at(&cube, Point3::new(1.0, 1.0, 1.0));
    for radii in [
        [0.2, 0.2, 0.2],
        [0.1, 0.15, 0.2],
        [0.05, 0.2, 0.3],
        [0.1, 0.1, 0.25],
    ] {
        let specs = edges_at(&cube, c)
            .into_iter()
            .zip(radii)
            .map(|(edge, r)| FilletSpec {
                edge,
                radius: FilletRadius::Constant(r),
            })
            .collect();
        run(
            &format!("cube_corner_{radii:?}"),
            &cube,
            specs,
            ContinuityLevel::G2,
        );
    }
    // All twelve cube edges.
    let specs = cube
        .edges
        .keys()
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    run("cube_12_edges", &cube, specs, ContinuityLevel::G2);
    // Notched saddle.
    let notch = Transform3::translation(0.6, 0.6, 0.6).compose(&Transform3::scale_uniform(0.6));
    let tool = BRep::unit_cube(tol)
        .expect("cube")
        .transform(&notch)
        .expect("tf");
    let body = popkern_bool::subtraction(BRep::unit_cube(tol).expect("cube"), tool, tol, None)
        .expect("notch");
    let s = vertex_at(&body, Point3::new(0.6, 0.6, 0.6));
    for (label, radii, cont) in [
        ("saddle_g1", [0.06, 0.08, 0.10], ContinuityLevel::G1),
        ("saddle_g2", [0.06, 0.08, 0.10], ContinuityLevel::G2),
        ("saddle_uniform_g2", [0.08, 0.08, 0.08], ContinuityLevel::G2),
    ] {
        let specs = edges_at(&body, s)
            .into_iter()
            .zip(radii)
            .map(|(edge, r)| FilletSpec {
                edge,
                radius: FilletRadius::Constant(r),
            })
            .collect();
        run(label, &body, specs, cont);
    }
}

#[test]
#[ignore = "diagnostic probe"]
fn census_pyramid_all_edges() {
    let tol = Tolerance::default();
    for &n in &[3usize, 4, 5] {
        for &h in &[0.8, 1.2] {
            let body = BRep::n_gonal_pyramid(n, 1.0, h, tol).expect("pyramid");
            let specs = body
                .edges
                .keys()
                .map(|edge| FilletSpec {
                    edge,
                    radius: FilletRadius::Constant(0.1),
                })
                .collect();
            run(
                &format!("pyr{n}_h{h}_all"),
                &body,
                specs,
                ContinuityLevel::G2,
            );
        }
    }
}

#[test]
#[ignore = "diagnostic probe"]
fn census_g1_requests() {
    pyramid("pyr3_h1.2_r0.16_g1", 3, 1.2, &[0.16], ContinuityLevel::G1);
    pyramid("pyr4_h1.2_r0.16_g1", 4, 1.2, &[0.16], ContinuityLevel::G1);
    pyramid("pyr5_h1.2_r0.16_g1", 5, 1.2, &[0.16], ContinuityLevel::G1);
    pyramid(
        "pyr4_h0.8_mixed_g1",
        4,
        0.8,
        &[0.04, 0.06, 0.05, 0.07],
        ContinuityLevel::G1,
    );
}
