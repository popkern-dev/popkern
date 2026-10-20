//! Curvature maps with a stated range (PK-504).

use super::*;
use popkern_core::tolerance::Tolerance;

fn map(brep: &BRep, kind: CurvatureKind, range: CurvatureRange) -> (Vec<u8>, CurvatureMapReport) {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::display())
        .expect("mesh");
    let mut camera = Camera::default_iso();
    camera.fit_to_bounds(glam::Vec3::splat(-1.0), glam::Vec3::splat(1.0));
    render_curvature_map(
        brep,
        &mesh,
        &camera,
        200,
        160,
        CurvatureMapOptions {
            kind,
            range,
            legend: false,
        },
    )
}

#[test]
fn a_spheres_map_reads_its_curvature_on_every_pixel_and_clips_nothing() {
    let r = 0.8;
    let sphere = BRep::sphere(r, Tolerance::default()).expect("sphere");
    let (_, report) = map(&sphere, CurvatureKind::Gaussian, CurvatureRange::default());
    let k = 1.0 / (r * r);
    assert!(report.pixels > 2000, "{report:?}");
    assert!(report.unevaluated * 100 < report.pixels, "{report:?}");
    assert!((report.sampled_min - k).abs() < 1e-6 * k, "{report:?}");
    assert!((report.sampled_max - k).abs() < 1e-6 * k, "{report:?}");
    assert!(!report.fixed && report.clipped == 0);
    assert!((report.half_range - report.sampled_max).abs() < 1e-12);
    // A fixed range under the value clips every evaluated pixel, and
    // says so.
    let (_, fixed) = map(&sphere, CurvatureKind::Gaussian, CurvatureRange::Fixed(1.0));
    assert!(fixed.fixed);
    assert_eq!(fixed.clipped, fixed.pixels - fixed.unevaluated);
}

#[test]
fn mean_curvature_is_signed_against_the_outward_normal() {
    // Convex surfaces bend away from their outward normal: H < 0 on a
    // sphere (-1/r) and on a cylinder's wall (-1/(2r)), whatever way
    // their charts face.
    let sphere = BRep::sphere(0.8, Tolerance::default()).expect("sphere");
    let (_, s) = map(&sphere, CurvatureKind::Mean, CurvatureRange::default());
    assert!((s.sampled_max + 1.0 / 0.8).abs() < 1e-6, "{s:?}");
    assert!((s.sampled_min + 1.0 / 0.8).abs() < 1e-6, "{s:?}");
    let cylinder = BRep::cylinder(0.5, 1.0, Tolerance::default()).expect("cylinder");
    let (_, c) = map(&cylinder, CurvatureKind::Mean, CurvatureRange::default());
    // Caps are flat (0); the wall reads -1/(2r) = -1.
    assert!((c.sampled_min + 1.0).abs() < 1e-6, "{c:?}");
    assert!(c.sampled_max.abs() < 1e-9, "{c:?}");
}

#[test]
fn the_legend_draws_the_bar_from_blue_to_red() {
    let sphere = BRep::sphere(0.8, Tolerance::default()).expect("sphere");
    let mesh = popkern_mesh::tessellate(&sphere, &popkern_mesh::TessellationOptions::display())
        .expect("mesh");
    let mut camera = Camera::default_iso();
    camera.fit_to_bounds(glam::Vec3::splat(-1.0), glam::Vec3::splat(1.0));
    let (w, h) = (240u32, 200u32);
    let (rgba, _) = render_curvature_map(
        &sphere,
        &mesh,
        &camera,
        w,
        h,
        CurvatureMapOptions {
            kind: CurvatureKind::Gaussian,
            range: CurvatureRange::default(),
            legend: true,
        },
    );
    let at = |x: u32, y: u32| {
        let i = ((y * w + x) * 4) as usize;
        [rgba[i], rgba[i + 1], rgba[i + 2]]
    };
    let x = w - 6;
    let (top, bottom) = (at(x, 14), at(x, h - 14));
    assert!(top[0] > top[2] + 50, "red at the top: {top:?}");
    assert!(bottom[2] > bottom[0] + 50, "blue at the bottom: {bottom:?}");
}
