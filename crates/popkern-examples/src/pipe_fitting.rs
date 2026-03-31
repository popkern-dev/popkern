//! Pipe fitting: a pipe swept along a 90-degree bend path.
//!
//! Demonstrates the `pipe` operation from popkern-sweep.

mod helpers;

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_sweep::{PipeRadius, SweepOptions, pipe};

fn main() {
    println!("=== Pipe Fitting ===");
    let tol = Tolerance::default();

    // Build a quarter-circle path in the XZ plane (90-degree bend).
    // Using a degree-2 rational NURBS arc: (0,0,0) -> (1,0,0) via (0,0,1).
    // This is a standard quarter-circle: center at (1,0,1), radius 1,
    // but let's use: start (0,0,0), mid-control (0,0,1), end (1,0,1).
    // Actually, a proper 90-degree arc needs a rational Bezier:
    //   P0 = (0,0,0), P1 = (0,0,1), P2 = (1,0,1)
    //   w = [1, 1/sqrt(2), 1]
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let path = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
        ],
        weights: vec![1.0, w, 1.0],
    };

    let radius = &PipeRadius::Constant(0.3);
    let options = SweepOptions::default();

    match pipe(&path, radius, &options, tol, None) {
        Ok(result) => {
            let brep = &result.brep;
            println!(
                "  Pipe: {} faces, {} edges, {} vertices",
                brep.faces.len(),
                brep.edges.len(),
                brep.vertices.len()
            );
            if !result.warnings.is_empty() {
                println!("  Warnings: {:?}", result.warnings);
            }
            helpers::export(brep, "Pipe Fitting", "pipe_fitting.step");
            helpers::try_render(brep, "pipe_fitting.png");
            helpers::render_multi_angle(brep, "pipe_fitting");
        }
        Err(e) => {
            eprintln!("  Pipe operation failed: {e}");
        }
    }

    println!("=== Done ===");
}
