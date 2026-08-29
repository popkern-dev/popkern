//! Offset/shell perf measurement harness (PK-329): per-fixture,
//! per-stage wall-clock for the pipeline over the corpus. Run in
//! release on an idle machine:
//!
//! ```bash
//! cargo run --release --bin offset_perf_probe
//! ```
//!
//! Times, per corpus body: `thicken` at the campaign-gate t = 0.02 with
//! per-stage durations from the progress tracker plus a standalone
//! re-run of the trim-stage SI scan and serve-gate validation for
//! attribution, `max_safe_offset` in both directions, the plane-only
//! cube reference, and the §8 plane-dominant ratio. These are the
//! numbers the perf-battery baselines in
//! `tests/perf_budget.rs` cite — re-run this binary when re-baselining
//! a budget row (the PK-190 discipline: a budget is `measured × 5`,
//! and every re-baseline commit names its measured source).

use popkern_core::progress::{CaptureReporter, OperationControl};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::corpus;
use popkern_shell::{ShellDirection, ShellOptions, max_safe_offset, thicken_with_control};
use std::time::Instant;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn main() {
    let fixtures = corpus(tol()).expect("corpus builds");

    // Plane-only reference: unit-cube shell (top face removed), the
    // "milliseconds on the cube battery" baseline from design §8.
    let cube = BRep::unit_cube(tol()).expect("cube");
    let top = cube
        .faces
        .keys()
        .max_by(|a, b| {
            let ca = cube.face_outer_loop_positions(*a);
            let cb = cube.face_outer_loop_positions(*b);
            let za = ca.iter().map(|p| p.z).sum::<f64>() / ca.len().max(1) as f64;
            let zb = cb.iter().map(|p| p.z).sum::<f64>() / cb.len().max(1) as f64;
            za.total_cmp(&zb)
        })
        .expect("cube has faces");
    let t0 = Instant::now();
    let r = popkern_shell::shell(
        &cube,
        0.05,
        ShellDirection::Inward,
        &[top],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    println!(
        "cube_shell_inward_0.05: {:.1} ms ({})",
        t0.elapsed().as_secs_f64() * 1e3,
        if r.is_ok() { "served" } else { "declined" }
    );

    // Plane-dominant ratio (design §8): cube thicken (pure plane
    // pipeline) vs l_bracket thicken (plane-dominant with concave
    // joins), best-of-5.
    let best_of = |n: usize, mut op: Box<dyn FnMut()>| -> f64 {
        op();
        (0..n)
            .map(|_| {
                let t0 = Instant::now();
                op();
                t0.elapsed().as_secs_f64() * 1e3
            })
            .fold(f64::MAX, f64::min)
    };
    let cube2 = BRep::unit_cube(tol()).expect("cube");
    let cube_ms = best_of(
        5,
        Box::new(move || {
            let _ = popkern_shell::thicken::thicken(&cube2, 0.02, &ShellOptions::default(), tol());
        }),
    );
    let lb = fixtures
        .iter()
        .find(|f| f.name == "l_bracket")
        .expect("l_bracket fixture");
    let lb_brep = lb.brep.clone();
    let lb_ms = best_of(
        5,
        Box::new(move || {
            let _ =
                popkern_shell::thicken::thicken(&lb_brep, 0.02, &ShellOptions::default(), tol());
        }),
    );
    println!(
        "plane ratio: cube thicken {cube_ms:.2} ms, l_bracket thicken {lb_ms:.2} ms, ratio {:.2}",
        lb_ms / cube_ms
    );

    for f in &fixtures {
        // max_safe_offset both directions.
        for dir in [ShellDirection::Inward, ShellDirection::Outward] {
            let t0 = Instant::now();
            let rep = max_safe_offset(&f.brep, dir, tol());
            let ms = t0.elapsed().as_secs_f64() * 1e3;
            match rep {
                Ok(rep) => println!(
                    "{}: max_safe {:?} = {:.4} in {:.1} ms",
                    f.name, dir, rep.max_distance, ms
                ),
                Err(e) => println!("{}: max_safe {:?} ERR {:.1} ms: {e}", f.name, dir, ms),
            }
        }

        // Thicken at the campaign-gate thickness with stage capture.
        let mut cap = CaptureReporter::default();
        let t0 = Instant::now();
        let result = thicken_with_control(
            &f.brep,
            0.02,
            &ShellOptions::default(),
            tol(),
            OperationControl::with_reporter(&mut cap),
        );
        let total = t0.elapsed().as_secs_f64() * 1e3;
        let verdict = match &result {
            Ok(_) => "served".to_string(),
            Err(e) => format!("declined: {e}"),
        };
        println!(
            "{}: thicken t=0.02 TOTAL {:.1} ms ({verdict})",
            f.name, total
        );
        let tracked: f64 = cap
            .events
            .iter()
            .filter_map(|e| e.duration)
            .map(|d| d.as_secs_f64() * 1e3)
            .sum();
        println!("    [tracked {tracked:.1} ms of {total:.1} ms]");
        for e in &cap.events {
            let ms = e.duration.map_or(f64::NAN, |d| d.as_secs_f64() * 1e3);
            if ms > 0.2 {
                println!("    {}: {ms:.1} ms", e.stage);
            }
        }
        // Attribute the untracked tail: the trim-stage SI scan and the
        // serve-gate validation, re-run standalone on the served body.
        if let Ok(res) = &result {
            let face_keys: Vec<_> = res.brep.faces.keys().collect();
            let t0 = Instant::now();
            let pairs = res
                .brep
                .collect_self_intersecting_pairs_among(&face_keys, true, true);
            let collect_ms = t0.elapsed().as_secs_f64() * 1e3;
            println!(
                "    [attrib] pair collection {collect_ms:.1} ms ({} pairs)",
                pairs.len()
            );
            let t0 = Instant::now();
            let scan =
                popkern_bool::ops::self_intersect::find_self_intersections(&res.brep, &tol());
            let si_ms = t0.elapsed().as_secs_f64() * 1e3;
            let t0 = Instant::now();
            let errors = res.brep.validate_collect_structural();
            let val_ms = t0.elapsed().as_secs_f64() * 1e3;
            println!(
                "    [attrib] SI scan {si_ms:.1} ms ({} curves), validate {val_ms:.1} ms ({} errs)",
                scan.map(|s| s.curves.len()).unwrap_or(0),
                errors.len()
            );
        }
    }
}
