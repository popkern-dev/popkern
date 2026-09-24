//! Declared-band derivation of `march_intersection` (PK-528): a
//! sub-precision endpoint residual is accepted only on fitted geometry
//! (a recorded band above `pos_tol`), and the curve declares 1.5× it,
//! capped under the Precision ceiling.

use super::*;
use popkern_core::geometry::Vec3;

const POS_TOL: f64 = 1e-5;

fn cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    }
}

fn floor_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

/// March the cylinder ∩ floor quarter circle from a start `offset`
/// off the cylinder (radially) to (0, 1, 0).
fn march(offset: f64, recorded: f64) -> Result<MarchedIntersection, SeamReason> {
    march_intersection(
        &cylinder(),
        &floor_plane(),
        Point3::new(1.0 + offset, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        POS_TOL,
        recorded,
    )
}

#[test]
fn exact_endpoints_keep_the_recorded_band() {
    let seam = march(0.0, POS_TOL).expect("exact quarter circle marches");
    assert_eq!(seam.band, POS_TOL);
    for p in &seam.samples[1..seam.samples.len() - 1] {
        let radial = (p.x * p.x + p.y * p.y).sqrt();
        assert!((radial - 1.0).abs() <= POS_TOL && p.z.abs() <= POS_TOL);
    }
}

#[test]
fn an_off_surface_endpoint_on_exact_geometry_declines() {
    assert_eq!(march(6e-5, POS_TOL).err(), Some(SeamReason::AnchorsOffSeam));
}

#[test]
fn a_residual_inside_the_recorded_band_keeps_the_band() {
    let seam = march(6e-5, 6.8e-5).expect("inside the recorded band");
    assert_eq!(seam.band, 6.8e-5);
}

#[test]
fn a_sub_precision_residual_on_fitted_geometry_declares_one_and_a_half_times_it() {
    let seam = march(6e-5, 5e-5).expect("sub-precision residual on fitted geometry");
    assert!((seam.band - 9e-5).abs() < 1e-12, "band {}", seam.band);
}

#[test]
fn a_residual_above_the_precision_cap_declines() {
    assert_eq!(march(2e-4, 5e-5).err(), Some(SeamReason::AnchorsOffSeam));
}
