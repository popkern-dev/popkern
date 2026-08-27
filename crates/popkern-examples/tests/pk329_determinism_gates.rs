//! PK-329 row 484: the offset/shell pipeline's parallel dispatch is
//! bit-identical at any thread count, pinned on the corpus.
//!
//! Every parallel site in the pipeline (per-face offset surfaces,
//! per-vertex and per-edge recompute, join classification, collapse
//! detection, the tier-3 fold scan, and the trim-stage SI pair scan)
//! follows the PK-277 ordered-collect discipline: a deterministic task
//! list, pure per-task reads, positional collection. These pins hold
//! the whole pipeline to that contract — thread pools of 1, 2, and 8
//! must produce byte-identical bodies and warnings, and a repeat run
//! must reproduce itself.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::{
    bracket_bore_countersink, filleted_cube, sphere_capped_rod, wavy_plate,
};
use popkern_shell::thicken::thicken;
use popkern_shell::{ShellOptions, ShellResult};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Structural + geometric bitwise comparison in handle-key order
/// (template: `boolean_determinism::assert_breps_identical`, PK-300).
fn assert_breps_identical(lhs: &BRep, rhs: &BRep, label: &str) {
    let keys = |b: &BRep| {
        (
            b.vertices
                .keys()
                .map(|k| format!("{k:?}"))
                .collect::<Vec<_>>(),
            b.edges.keys().map(|k| format!("{k:?}")).collect::<Vec<_>>(),
            b.faces.keys().map(|k| format!("{k:?}")).collect::<Vec<_>>(),
            b.loops.keys().map(|k| format!("{k:?}")).collect::<Vec<_>>(),
        )
    };
    assert_eq!(keys(lhs), keys(rhs), "{label}: handle sequences differ");

    for (lv, rv) in lhs.vertices.values().zip(rhs.vertices.values()) {
        let bits =
            |p: &popkern_core::geometry::Point3| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
        assert_eq!(
            bits(&lv.position),
            bits(&rv.position),
            "{label}: vertex positions differ"
        );
    }
    for ((le_h, le), (_, re)) in lhs.edges.iter().zip(rhs.edges.iter()) {
        match (&le.curve, &re.curve) {
            (None, None) => {}
            (Some(lc), Some(rc)) => {
                assert_eq!(
                    (lc.t_start.to_bits(), lc.t_end.to_bits()),
                    (rc.t_start.to_bits(), rc.t_end.to_bits()),
                    "{label}: edge {le_h:?} parameter window differs"
                );
                for k in 0..=8 {
                    let t = lc.t_start + (lc.t_end - lc.t_start) * (k as f64) / 8.0;
                    let (lp, rp) = (lc.geometry.point_at(t), rc.geometry.point_at(t));
                    assert_eq!(
                        (lp.x.to_bits(), lp.y.to_bits(), lp.z.to_bits()),
                        (rp.x.to_bits(), rp.y.to_bits(), rp.z.to_bits()),
                        "{label}: edge {le_h:?} curve sample {k} differs"
                    );
                }
            }
            _ => panic!("{label}: edge {le_h:?} curve presence differs"),
        }
    }
    for ((lf_h, lf), (_, rf)) in lhs.faces.iter().zip(rhs.faces.iter()) {
        match (&lf.surface, &rf.surface) {
            (None, None) => {}
            (Some(ls), Some(rs)) => {
                assert_eq!(
                    ls.orientation, rs.orientation,
                    "{label}: face {lf_h:?} orientation flag differs"
                );
                let ((lu0, lu1), (lv0, lv1)) = ls.geometry.parameter_range();
                for (i, j) in [(1, 1), (1, 3), (3, 1), (2, 2)] {
                    let sample = |s: &popkern_core::geometry::surface::SurfaceGeometry,
                                  (u0, u1): (f64, f64),
                                  (v0, v1): (f64, f64)| {
                        let cl = |lo: f64, hi: f64, f: f64| {
                            if lo.is_finite() && hi.is_finite() {
                                lo + (hi - lo) * f
                            } else {
                                f
                            }
                        };
                        s.point_at(cl(u0, u1, i as f64 / 4.0), cl(v0, v1, j as f64 / 4.0))
                    };
                    let lp = sample(&ls.geometry, (lu0, lu1), (lv0, lv1));
                    let rp = sample(&rs.geometry, (lu0, lu1), (lv0, lv1));
                    assert_eq!(
                        (lp.x.to_bits(), lp.y.to_bits(), lp.z.to_bits()),
                        (rp.x.to_bits(), rp.y.to_bits(), rp.z.to_bits()),
                        "{label}: face {lf_h:?} surface sample ({i},{j}) differs"
                    );
                }
            }
            _ => panic!("{label}: face {lf_h:?} surface presence differs"),
        }
    }
}

fn assert_results_identical(lhs: &ShellResult, rhs: &ShellResult, label: &str) {
    assert_breps_identical(&lhs.brep, &rhs.brep, label);
    assert_eq!(
        format!("{:?}", lhs.warnings),
        format!("{:?}", rhs.warnings),
        "{label}: warnings differ"
    );
}

fn thicken_in_pool(brep: &BRep, threads: usize) -> ShellResult {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("pool builds");
    pool.install(|| thicken(brep, 0.02, &ShellOptions::default(), tol()).expect("thicken serves"))
}

fn pin_body(brep: &BRep, label: &str) {
    let one = thicken_in_pool(brep, 1);
    let two = thicken_in_pool(brep, 2);
    let eight = thicken_in_pool(brep, 8);
    assert_results_identical(&one, &two, &format!("{label}: 1 vs 2 threads"));
    assert_results_identical(&one, &eight, &format!("{label}: 1 vs 8 threads"));
    let repeat = thicken_in_pool(brep, 8);
    assert_results_identical(&eight, &repeat, &format!("{label}: repeat at 8 threads"));
}

#[test]
fn thread_count_does_not_change_filleted_cube_thicken() {
    // Spheres, cylinders, planes: per-face offsets, join
    // classification, and the analytic-certificate pair scan.
    let f = filleted_cube(tol()).expect("fixture builds");
    pin_body(&f.brep, "filleted_cube");
}

#[test]
fn thread_count_does_not_change_wavy_plate_thicken() {
    // Free-form NURBS: fitted offsets, the parallel general-SSI edge
    // recompute (the banded fit), and the tier-3 fold scan.
    let f = wavy_plate(tol()).expect("fixture builds");
    pin_body(&f.brep, "wavy_plate");
}

#[test]
fn thread_count_does_not_change_bracket_thicken() {
    // Bore + countersink: cones (windowed nappe certificates) and the
    // coaxial-profile edge arms.
    let f = bracket_bore_countersink(tol()).expect("fixture builds");
    pin_body(&f.brep, "bracket_bore_countersink");
}

#[test]
fn thread_count_does_not_change_rod_thicken() {
    // Boolean-built rod: concentric sphere/cylinder certificate class.
    let f = sphere_capped_rod(tol()).expect("fixture builds");
    pin_body(&f.brep, "sphere_capped_rod");
}
