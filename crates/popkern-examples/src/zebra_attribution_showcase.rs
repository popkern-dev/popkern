//! Zebra stripe attribution (PK-503): which stripe defects on the
//! setback vertex caps belong to the SURFACE and which to the DISPLAY
//! MESH.
//!
//! Each fixture is rendered once into a deferred zebra frame and read
//! in both normal sources — interpolated display-mesh normals (what
//! the playground shows) and the faces' own normals at the pixel's ray
//! hit. Per seam class (cap–cap, cap rim = cap–strip, strip–support,
//! strip–strip) the image-space seam readings are reported for both
//! views: the NORMAL jump (radians, both sides extrapolated to the
//! seam) is the G1 break, the reflection coordinate's SLOPE jump (per
//! pixel) the G2 break. A reading present in the surface
//! view is the geometry's; one present only in the mesh view is the
//! tessellation's.
//!
//! Outputs: `examples/output/zebra_attribution_<fixture>_{mesh,surface}.png`
//! and the measurement table on stdout.

use std::path::Path;

use popkern_examples::zebra_attribution::{attribute, attribution_fixtures};
use popkern_viz::zebra::ZebraNormals;

fn main() {
    println!("=== Zebra attribution (PK-503) ===");
    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/output");
    std::fs::create_dir_all(&out_dir).expect("output dir");
    let (size, stripes) = (640u32, 16u32);
    for fx in attribution_fixtures() {
        let at = attribute(&fx, size, stripes);
        for (tag, normals) in [
            ("mesh", ZebraNormals::Mesh),
            ("surface", ZebraNormals::Surface),
        ] {
            let path = out_dir.join(format!("zebra_attribution_{}_{tag}.png", fx.name));
            popkern_viz::headless::write_png(&path, size, size, &at.frame.shade(normals, stripes))
                .expect("png");
        }
        let report = &at.disagreement;
        println!(
            "\n{}: {} faces, {} cap faces; normals mesh-vs-surface max {:.3}° mean {:.4}° p99 {:.3}°, {} px change band ({} px, {} ray hits)",
            fx.name,
            fx.brep.faces.len(),
            fx.caps.len(),
            report.max_angle_rad.to_degrees(),
            report.mean_angle_rad.to_degrees(),
            report.p99_angle_rad.to_degrees(),
            report.band_mismatch_pixels,
            report.evaluated_pixels,
            report.ray_hit_pixels,
        );
        println!(
            "  {:<14} {:>5} {:>3} | {:>9} {:>9} | {:>9} {:>9} | {:>9} {:>9} | {:>9} {:>9}",
            "seam class",
            "rows",
            "sm",
            "surf nj50",
            "mesh nj50",
            "surf nj95",
            "mesh nj95",
            "surf sl50",
            "mesh sl50",
            "surf sl95",
            "mesh sl95"
        );
        for (class, s, m) in &at.classes {
            println!(
                "  {:<14} {:>5} {:>3} | {:>9.2e} {:>9.2e} | {:>9.2e} {:>9.2e} | {:>9.2e} {:>9.2e} | {:>9.2e} {:>9.2e}",
                class.label(),
                s.rows,
                s.seams,
                s.normal_p50,
                m.normal_p50,
                s.normal_p95,
                m.normal_p95,
                s.slope_p50,
                m.slope_p50,
                s.slope_p95,
                m.slope_p95
            );
        }
    }
}
