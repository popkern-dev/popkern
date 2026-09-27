//! Tests for the attributable zebra frame: on analytic bodies the
//! surface view reproduces the exact normal and is independent of the
//! tessellation, while the mesh view carries the facets' error.

use super::locate::unwrap_near;
use super::*;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_mesh::TessellationOptions;

fn coarse() -> TessellationOptions {
    TessellationOptions {
        chord_height: 0.08,
        min_subdivisions: 2,
        max_subdivisions: 6,
        max_facet_angle: std::f64::consts::PI,
    }
}

fn fine() -> TessellationOptions {
    TessellationOptions {
        chord_height: 5e-4,
        min_subdivisions: 8,
        max_subdivisions: 128,
        max_facet_angle: 3f64.to_radians(),
    }
}

fn camera() -> Camera {
    Camera {
        eye: glam::Vec3::new(2.2, -1.8, 1.4),
        target: glam::Vec3::ZERO,
        up: glam::Vec3::Z,
        fov_y: std::f32::consts::FRAC_PI_4,
        aspect: 4.0 / 3.0,
        near: 0.01,
        far: 100.0,
    }
}

/// The largest angle between the frame's surface normals and `truth`
/// evaluated at each sample's surface point.
fn worst_surface_error(
    brep: &BRep,
    frame: &ZebraFrame,
    truth: impl Fn([f64; 3]) -> [f64; 3],
) -> (f64, usize) {
    let mut worst = 0.0f64;
    let mut n = 0usize;
    for s in frame.samples.iter().flatten() {
        let Some(sn) = s.surface_normal else { continue };
        let Some(sof) = brep.faces.get(s.face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let p = sof.geometry.point_at(s.uv.0, s.uv.1);
        let t = truth([p.x, p.y, p.z]);
        worst = worst.max(angle_between(sn, t));
        n += 1;
    }
    (worst, n)
}

/// Pixels covered in both frames whose stripe band differs, as a
/// fraction of the pixels covered in both.
fn band_difference(
    a: &ZebraFrame,
    b: &ZebraFrame,
    normals: ZebraNormals,
    stripe_count: u32,
) -> f64 {
    let mut both = 0usize;
    let mut differ = 0usize;
    for (sa, sb) in a.samples.iter().zip(&b.samples) {
        if let (Some(sa), Some(sb)) = (sa, sb) {
            both += 1;
            if a.band(sa, normals, stripe_count) != b.band(sb, normals, stripe_count) {
                differ += 1;
            }
        }
    }
    differ as f64 / both.max(1) as f64
}

#[test]
fn sphere_surface_view_is_the_exact_normal_and_tessellation_independent() {
    let tol = Tolerance::default();
    let brep = BRep::sphere(0.8, tol).expect("sphere");
    let cam = camera();
    let coarse_mesh = popkern_mesh::tessellate(&brep, &coarse()).expect("coarse");
    let fine_mesh = popkern_mesh::tessellate(&brep, &fine()).expect("fine");
    let cf = rasterize_zebra_frame(&brep, &coarse_mesh, &cam, 240, 180);
    let ff = rasterize_zebra_frame(&brep, &fine_mesh, &cam, 240, 180);
    assert!(cf.mixed_facing_faces.is_empty());

    // Surface view == outward radial direction at the sample's surface
    // point, at any tessellation (the pole limit included).
    let radial = |p: [f64; 3]| unit(p).unwrap_or([0.0, 0.0, 1.0]);
    let (err, n) = worst_surface_error(&brep, &cf, radial);
    assert!(n > 5000, "sphere must cover the frame ({n} px)");
    assert!(err < 1e-9, "surface normal must be exact, worst {err:.3e}");

    // Every pixel's surface point is the ray hit (convex body: the ray
    // crosses the surface before reaching the inscribed facet).
    let report = cf.mode_disagreement(12);
    assert_eq!(report.evaluated_pixels, report.covered_pixels);
    assert_eq!(report.ray_hit_pixels, report.evaluated_pixels);
    let worst = report.worst.expect("worst pixel");
    assert!((worst.angle_rad - report.max_angle_rad).abs() < 1e-15);
    let fine_report = ff.mode_disagreement(12);
    assert!(
        fine_report.max_angle_rad < 0.5 * report.max_angle_rad,
        "refinement must shrink the disagreement ({:.3e} vs {:.3e})",
        fine_report.max_angle_rad,
        report.max_angle_rad
    );

    // Stripes in the surface view hardly move when the mesh changes;
    // in the mesh view they do. (The residue is silhouette and
    // chord-offset view vectors: both frames shade the FACET point.)
    let surface_shift = band_difference(&cf, &ff, ZebraNormals::Surface, 12);
    let mesh_shift = band_difference(&cf, &ff, ZebraNormals::Mesh, 12);
    eprintln!("sphere band shift coarse->fine: surface {surface_shift:.4}, mesh {mesh_shift:.4}");
    assert!(
        surface_shift < 1e-3 && mesh_shift > 0.1,
        "surface view must be tessellation-stable: {surface_shift:.4} vs mesh {mesh_shift:.4}"
    );
}

#[test]
fn cylinder_surface_view_is_exact_across_the_seam() {
    let tol = Tolerance::default();
    let brep = BRep::cylinder(0.5, 1.0, tol).expect("cylinder");
    let cam = Camera {
        target: glam::Vec3::new(0.0, 0.0, 0.5),
        ..camera()
    };
    let mesh = popkern_mesh::tessellate(&brep, &coarse()).expect("coarse");
    let frame = rasterize_zebra_frame(&brep, &mesh, &cam, 240, 180);
    // Lateral samples: outward horizontal radial; caps: ±z.
    let truth = |p: [f64; 3]| {
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        if r < 0.5 - 1e-6 {
            [0.0, 0.0, if p[2] > 0.5 { 1.0 } else { -1.0 }]
        } else {
            [p[0] / r, p[1] / r, 0.0]
        }
    };
    let (err, n) = worst_surface_error(&brep, &frame, truth);
    assert!(n > 3000, "cylinder must cover the frame ({n} px)");
    assert!(err < 1e-9, "surface normal must be exact, worst {err:.3e}");
    let report = frame.mode_disagreement(12);
    assert_eq!(report.ray_hit_pixels, report.evaluated_pixels);
    // Planes carry no facet error: the per-face breakdown attributes
    // the whole disagreement to the lateral wall.
    let lateral = report.per_face.first().expect("per-face");
    for other in &report.per_face[1..] {
        assert!(
            other.max_angle_rad < 1e-6,
            "cap faces must agree exactly, got {:.3e}",
            other.max_angle_rad
        );
    }
    assert!(lateral.max_angle_rad > 1e-3);
}

/// A torus's normal is not a linear function of position, so the
/// normalized-lerp mesh normal errs inside coarse facets (a sphere or
/// cylinder hides this: lerped radial normals ARE the radial direction
/// of the chord point). The comparison must attribute that error to
/// the mesh, and it must shrink under refinement.
#[test]
fn torus_mesh_view_error_is_attributed_to_the_tessellation() {
    let tol = Tolerance::default();
    let brep = BRep::torus(0.7, 0.25, tol).expect("torus");
    let cam = camera();
    let coarse_mesh = popkern_mesh::tessellate(&brep, &coarse()).expect("coarse");
    let fine_mesh = popkern_mesh::tessellate(&brep, &fine()).expect("fine");
    let cf = rasterize_zebra_frame(&brep, &coarse_mesh, &cam, 240, 180);
    let ff = rasterize_zebra_frame(&brep, &fine_mesh, &cam, 240, 180);
    let rc = cf.mode_disagreement(12);
    let rf = ff.mode_disagreement(12);
    eprintln!(
        "torus coarse max {:.3e} mean {:.3e} ray {} cp {} interp {}; fine max {:.3e} mean {:.3e}",
        rc.max_angle_rad,
        rc.mean_angle_rad,
        rc.ray_hit_pixels,
        rc.closest_point_pixels,
        rc.interpolated_pixels,
        rf.max_angle_rad,
        rf.mean_angle_rad
    );
    assert!(rc.mean_angle_rad > 4.0 * rf.mean_angle_rad);
    assert!(rc.band_mismatch_pixels > rf.band_mismatch_pixels);
    let surface_shift = band_difference(&cf, &ff, ZebraNormals::Surface, 12);
    let mesh_shift = band_difference(&cf, &ff, ZebraNormals::Mesh, 12);
    eprintln!("torus band shift coarse->fine: surface {surface_shift:.4}, mesh {mesh_shift:.4}");
    assert!(surface_shift < 2e-3 && mesh_shift > 0.1);
}

#[test]
fn periodic_uv_unwraps_onto_one_sheet() {
    let tau = std::f64::consts::TAU;
    assert!((unwrap_near(tau - 0.01, 0.02, Some(tau)) - (-0.01)).abs() < 1e-12);
    assert!((unwrap_near(0.01, tau - 0.02, Some(tau)) - (tau + 0.01)).abs() < 1e-12);
    assert_eq!(unwrap_near(5.0, 0.0, None), 5.0);
}

#[test]
fn shading_paints_both_bands_and_background() {
    let tol = Tolerance::default();
    let brep = BRep::sphere(0.8, tol).expect("sphere");
    let mesh = popkern_mesh::tessellate(&brep, &coarse()).expect("mesh");
    let frame = rasterize_zebra_frame(&brep, &mesh, &camera(), 160, 120);
    for normals in [ZebraNormals::Mesh, ZebraNormals::Surface] {
        let img = frame.shade(normals, 8);
        assert_eq!(img.len(), 160 * 120 * 4);
        let dark = img.chunks(4).filter(|p| p[0] == DARK).count();
        let light = img.chunks(4).filter(|p| p[0] == LIGHT).count();
        let bg = img.chunks(4).filter(|p| *p == BACKGROUND).count();
        assert!(dark > 100 && light > 100 && bg > 100, "{normals:?}");
    }
    assert!(
        frame
            .reflection_coordinate(80, 60, ZebraNormals::Surface)
            .is_some()
    );
    assert!(
        frame
            .reflection_coordinate(0, 0, ZebraNormals::Surface)
            .is_none()
    );
    assert!(
        frame
            .reflection_coordinate(500, 0, ZebraNormals::Mesh)
            .is_none()
    );
}
