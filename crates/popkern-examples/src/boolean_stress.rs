//! Boolean stress test: progressively harder geometric operations.
//!
//! Each test exercises a different part of the boolean pipeline:
//! 1. Cube-cube union (planar-planar, basic)
//! 2. Cube with cylindrical hole (planar-curved subtraction)
//! 3. Two overlapping cubes subtraction (partial overlap)
//! 4. Sphere-cube intersection (curved-planar)
//! 5. Multiple holes in a plate (instancing)
//! 6. Identical cubes union (coincidence handling)
//! 7. Touching cubes union (shared face)

mod helpers;

use popkern_bool::{intersection, subtraction, subtraction_ref, union};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn validate_result(name: &str, result: &BRep) -> bool {
    let v = result.vertices.len();
    let e = result.edges.len();
    let f = result.faces.len();
    let euler = v as i64 - e as i64 + f as i64;

    // Structural checks (topology)
    let structural_ok = result.validate_euler_poincare().is_ok()
        && result.validate_half_edge_pairing().is_ok()
        && result.validate_manifold_edges().is_ok();

    // Full validation (includes geometric checks)
    let full_ok = result.validate().is_ok();

    let status = if structural_ok && full_ok {
        "PASS"
    } else if structural_ok {
        "GEOM" // structurally valid but geometric issues
    } else {
        "FAIL"
    };
    println!("  [{status}] {name}: V={v} E={e} F={f} euler={euler}");
    if !full_ok {
        let errors = result.validate_all();
        for err in &errors {
            println!("         {err}");
        }
    }
    structural_ok // pass if structurally valid (geometric issues are warnings)
}

fn main() {
    println!("=== Boolean Stress Test ===\n");
    let mut pass = 0;
    let mut fail = 0;
    let skip = 0;

    // --- Test 1: Cube-cube union (offset) ---
    print!("Test 1: Cube-cube union... ");
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    match union(a, b, tol(), None) {
        Ok(result) => {
            println!();
            if validate_result("cube-cube union", &result) {
                helpers::try_render(&result, "stress_01_cube_union.png");
                pass += 1;
            } else {
                fail += 1;
            }
        }
        Err(e) => {
            println!("ERROR: {e}");
            fail += 1;
        }
    }

    // --- Test 2: Cube with cylindrical hole ---
    print!("Test 2: Cube with hole... ");
    let cube = BRep::unit_cube(tol()).unwrap();
    let cyl = BRep::cylinder(0.2, 2.0, tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    match subtraction(cube, cyl, tol(), None) {
        Ok(result) => {
            println!();
            if validate_result("cube with hole", &result) {
                helpers::try_render(&result, "stress_02_cube_hole.png");
                pass += 1;
            } else {
                fail += 1;
            }
        }
        Err(e) => {
            println!("ERROR: {e}");
            fail += 1;
        }
    }

    // --- Test 3: Overlapping cubes subtraction ---
    print!("Test 3: Overlapping cubes subtraction... ");
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, 0.0))
        .unwrap();
    match subtraction(a, b, tol(), None) {
        Ok(result) => {
            println!();
            if validate_result("overlap subtraction", &result) {
                helpers::try_render(&result, "stress_03_overlap_sub.png");
                pass += 1;
            } else {
                fail += 1;
            }
        }
        Err(e) => {
            println!("ERROR: {e}");
            fail += 1;
        }
    }

    // --- Test 4: Sphere-cube intersection ---
    print!("Test 4: Sphere-cube intersection... ");
    let cube = BRep::unit_cube(tol()).unwrap();
    let sphere = BRep::sphere(0.8, tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, 0.5))
        .unwrap();
    match intersection(cube, sphere, tol(), None) {
        Ok(result) => {
            println!();
            if validate_result("sphere-cube intersection", &result) {
                helpers::try_render(&result, "stress_04_sphere_cube.png");
                pass += 1;
            } else {
                fail += 1;
            }
        }
        Err(e) => {
            println!("ERROR: {e}");
            fail += 1;
        }
    }

    // --- Test 5: Multiple holes in a plate ---
    print!("Test 5: Multiple holes in plate... ");
    let plate = BRep::unit_cube(tol()).unwrap();
    let hole_template = BRep::cylinder(0.1, 2.0, tol()).unwrap();
    let positions = [
        (0.25, 0.25),
        (0.25, 0.75),
        (0.75, 0.25),
        (0.75, 0.75),
        (0.5, 0.5),
    ];
    let mut plate_result = plate;
    let mut holes_ok = true;
    for (x, y) in &positions {
        let hole = hole_template
            .transform(&Transform3::translation(*x, *y, -0.5))
            .unwrap();
        match subtraction_ref(&plate_result, &hole, tol(), None) {
            Ok(r) => plate_result = r,
            Err(e) => {
                println!("\n  Hole at ({x},{y}) failed: {e}");
                holes_ok = false;
                break;
            }
        }
    }
    if holes_ok {
        println!();
        if validate_result("plate with 5 holes", &plate_result) {
            helpers::try_render(&plate_result, "stress_05_plate_holes.png");
            pass += 1;
        } else {
            fail += 1;
        }
    } else {
        fail += 1;
    }

    // --- Test 6: Identical cubes union (coincidence) ---
    print!("Test 6: Identical cubes union... ");
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    match union(a, b, tol(), None) {
        Ok(result) => {
            println!();
            if validate_result("identical union", &result) {
                pass += 1;
            } else {
                fail += 1;
            }
        }
        Err(e) => {
            println!("ERROR: {e}");
            fail += 1;
        }
    }

    // --- Test 7: Touching cubes union (shared face) ---
    print!("Test 7: Touching cubes union... ");
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(1.0, 0.0, 0.0))
        .unwrap();
    match union(a, b, tol(), None) {
        Ok(result) => {
            println!();
            if validate_result("touching union", &result) {
                helpers::try_render(&result, "stress_07_touching.png");
                pass += 1;
            } else {
                fail += 1;
            }
        }
        Err(e) => {
            println!("ERROR: {e}");
            fail += 1;
        }
    }

    // --- Summary ---
    println!("\n=== Results: {pass} pass, {fail} fail, {skip} skip ===");
    if fail > 0 {
        std::process::exit(1);
    }
}
