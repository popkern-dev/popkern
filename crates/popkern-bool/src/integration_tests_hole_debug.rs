//! Debug tests for cylinder-through-cube hole subtraction.
//!
//! The stress test revealed that cube-cylinder subtraction at tight tolerance
//! (1e-6) produces invalid topology, while at coarse tolerance (1e-2) it works.
//! These tests trace the issue.

use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn default_tol() -> Tolerance {
    Tolerance::default() // 1e-6
}

fn coarse_tol() -> Tolerance {
    Tolerance::default().with_linear(1e-2)
}

fn medium_tol() -> Tolerance {
    Tolerance::default().with_linear(1e-4)
}

fn make_cube_with_hole(tol: Tolerance) -> Result<BRep, crate::error::BoolError> {
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.2, 2.0, tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    crate::subtraction(cube, cyl, tol, None)
}

#[test]
fn hole_at_coarse_tol_works() {
    let result = make_cube_with_hole(coarse_tol()).expect("coarse tol should work");
    let v = result.vertices.len();
    let e = result.edges.len();
    let f = result.faces.len();
    eprintln!("[coarse] V={v} E={e} F={f}");
    assert!(f >= 6, "should have at least 6 faces");
    let inner_loops: usize = result.faces.values().map(|f| f.inner_loops.len()).sum();
    eprintln!("[coarse] inner_loops={inner_loops}");
}

#[test]
fn hole_at_medium_tol() {
    let result = make_cube_with_hole(medium_tol());
    match result {
        Ok(brep) => {
            let v = brep.vertices.len();
            let e = brep.edges.len();
            let f = brep.faces.len();
            let errors = brep.validate_all();
            eprintln!("[medium] V={v} E={e} F={f} errors={}", errors.len());
            for err in &errors {
                eprintln!("[medium]   {err}");
            }
        }
        Err(e) => {
            eprintln!("[medium] FAILED: {e}");
        }
    }
}

#[test]
fn hole_at_default_tol() {
    let result = make_cube_with_hole(default_tol());
    match result {
        Ok(brep) => {
            let v = brep.vertices.len();
            let e = brep.edges.len();
            let f = brep.faces.len();
            let errors = brep.validate_all();
            eprintln!("[default] V={v} E={e} F={f} errors={}", errors.len());
            for err in &errors {
                eprintln!("[default]   {err}");
            }

            // Check intersection curve count
            let tol = default_tol();
            let cube = BRep::unit_cube(tol).unwrap();
            let cyl = BRep::cylinder(0.2, 2.0, tol)
                .unwrap()
                .transform(&Transform3::translation(0.5, 0.5, -0.5))
                .unwrap();
            let isect = crate::intersect::find_intersections(&cube, &cyl, &tol);
            match isect {
                Ok(result) => {
                    eprintln!(
                        "[default] intersection: {} curves, {} coincident overlaps",
                        result.records.len(),
                        result.coincidence.face_overlaps.len()
                    );
                    for (i, rec) in result.records.iter().enumerate() {
                        let (t0, t1) = rec.curve.curve_3d.parameter_range();
                        use popkern_core::geometry::CurveEval;
                        let p0 = rec.curve.curve_3d.point_at(t0);
                        let p1 = rec.curve.curve_3d.point_at(t1);
                        let closed = rec.curve.is_closed;
                        eprintln!(
                            "[default]   curve {i}: face_a={:?} face_b={:?} closed={closed}",
                            rec.face_a, rec.face_b
                        );
                        eprintln!(
                            "[default]     start=({:.3},{:.3},{:.3}) end=({:.3},{:.3},{:.3})",
                            p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                        );
                    }
                }
                Err(diag) => {
                    eprintln!("[default] intersection FAILED: {:?}", diag.error);
                }
            }
        }
        Err(e) => {
            eprintln!("[default] FAILED: {e}");
        }
    }
}
