//! Boolean pipeline repeat-call determinism pin (PK-300).
//!
//! PK-277's tessellation identity pins caught two real pre-existing
//! hash-order nondeterminisms (Bowyer–Watson cavity refill, crack-
//! zipper candidate ties), both invisible to topology-count goldens.
//! This pin holds the boolean pipeline to the same standard: running
//! the same operation twice on the same operands must produce
//! bit-identical BReps — same entity counts, same handle sequences,
//! same vertex positions, same sampled edge curves and face surfaces.
//! A `HashMap` iterated anywhere in split/classify/merge/postprocess
//! in an order that feeds output construction or tie-breaking would
//! break this.

use popkern_core::geometry::{CurveEval, SurfaceEval, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// The PK-295 served-subtraction operands: a rim-filleted cylinder
/// (curved walls, torus strip, closed rims) and a strip-crossing
/// cuboid — exercises curved SSI, split, classification, merge, and
/// the postprocess passes.
fn operands() -> (BRep, BRep) {
    let cyl = BRep::cylinder(1.0, 2.0, tol()).expect("cylinder builds");
    let rim: Vec<_> = cyl
        .edges
        .iter()
        .filter(|(_, d)| {
            let he = &cyl.half_edges[d.half_edge];
            let a = cyl.vertices[he.origin].position;
            let b = cyl.vertices[cyl.half_edges[he.twin].origin].position;
            let f1 = cyl.loops[he.loop_].face;
            let f2 = cyl.loops[cyl.half_edges[he.twin].loop_].face;
            f1 != f2 && (a.z - 2.0).abs() < 1e-9 && (b.z - 2.0).abs() < 1e-9
        })
        .map(|(e, _)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.2),
        })
        .collect();
    let part = fillet_edges(&cyl, &rim, &FilletOptions::default(), tol(), None)
        .expect("rim fillet completes")
        .brep;
    let cutter = BRep::unit_cube(tol())
        .expect("cube builds")
        .transform(&Transform3::scale_uniform(1.8))
        .expect("scales")
        .transform(&Transform3::translation(-0.9, -0.9, 0.6))
        .expect("translates");
    (part, cutter)
}

/// Structural + geometric bitwise comparison in handle-key order.
/// Handle sequences are compared via their Debug form: identical
/// construction order yields identical slot/version sequences.
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

/// Repeat subtraction on identical operands must be bit-identical.
#[test]
fn repeat_subtraction_is_bit_identical() {
    let (part, cutter) = operands();
    let r1 = popkern_bool::subtraction_ref(&part, &cutter, tol(), None).expect("first subtraction");
    let r2 =
        popkern_bool::subtraction_ref(&part, &cutter, tol(), None).expect("second subtraction");
    assert_breps_identical(&r1, &r2, "repeat subtraction");
}

/// Repeat subtraction must also be thread-count independent end to
/// end (Stage 1 runs on rayon; later stages are sequential today —
/// this pins that contract for them too).
#[test]
fn thread_count_does_not_change_subtraction() {
    let (part, cutter) = operands();
    let single = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("1-thread pool builds")
        .install(|| {
            popkern_bool::subtraction_ref(&part, &cutter, tol(), None).expect("subtraction")
        });
    let many = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .expect("8-thread pool builds")
        .install(|| {
            popkern_bool::subtraction_ref(&part, &cutter, tol(), None).expect("subtraction")
        });
    assert_breps_identical(&single, &many, "1-thread vs 8-thread subtraction");
}
