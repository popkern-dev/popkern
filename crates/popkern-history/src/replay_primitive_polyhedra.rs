//! The flat-faced primitives (PK-1083): a regular prism (an extruded
//! polygon, as a box is), a regular pyramid and its frustum, and a
//! wedge, the last two built exactly by the core constructors from
//! their corners. All in the local frame: axis +z, corner 0 (or the
//! box corner) on +x.

use popkern_core::geometry::Point3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;
use popkern_sketch::profile::{AnalyticProfile, ProfileShape};

use super::{extruded, failed, invalid};
use crate::error::ReplayError;

/// The corners of a regular `sides`-gon at height `z`, corner 0 on
/// +x; `radius` to the corners, or to the flats when `flats`.
fn polygon(sides: u32, radius: f64, flats: bool, z: f64) -> Vec<Point3> {
    let n = sides as f64;
    let corner = if flats {
        radius / (std::f64::consts::PI / n).pcos()
    } else {
        radius
    };
    (0..sides)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n;
            Point3::new(corner * a.pcos(), corner * a.psin(), z)
        })
        .collect()
}

fn check_sides(shape: &str, sides: u32) -> Result<(), ReplayError> {
    if sides < 3 {
        return Err(invalid(shape, "sides", sides as f64, "needs at least 3"));
    }
    Ok(())
}

/// A regular prism standing on the origin.
pub(super) fn prism(
    shape: &str,
    sides: u32,
    radius: f64,
    height: f64,
    flats: bool,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    check_sides(shape, sides)?;
    let points = polygon(sides, radius, flats, 0.0)
        .iter()
        .map(|p| [p.x, p.y])
        .collect();
    let profile = AnalyticProfile::new(ProfileShape::Polygon { points });
    extruded(shape, &profile, 0.0, height, tolerance)
}

/// A regular pyramid (`r1` or `r2` 0: an apex) or its frustum.
pub(super) fn pyramid(
    shape: &str,
    sides: u32,
    (r1, r2): (f64, f64),
    height: f64,
    flats: bool,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    check_sides(shape, sides)?;
    let fail = |e: popkern_core::error::CoreError| failed(shape, "pyramid", e);
    match (r1 > 0.0, r2 > 0.0) {
        (true, true) => BRep::polygon_frustum(
            &polygon(sides, r1, flats, 0.0),
            &polygon(sides, r2, flats, height),
            tolerance,
        )
        .map_err(fail),
        (true, false) => BRep::n_gonal_pyramid(
            sides as usize,
            polygon_corner(sides, r1, flats),
            height,
            tolerance,
        )
        .map_err(fail),
        (false, true) => {
            // Apex down: the upright pyramid turned over about x (its
            // corner 0 stays on +x) and lifted by the height.
            let upright = BRep::n_gonal_pyramid(
                sides as usize,
                polygon_corner(sides, r2, flats),
                height,
                tolerance,
            )
            .map_err(fail)?;
            let over = Transform3::translation(0.0, 0.0, height)
                .compose(&Transform3::rotation_x(std::f64::consts::PI));
            upright
                .transform(&over)
                .map_err(|e| failed(shape, "pyramid", e))
        }
        (false, false) => Err(invalid(shape, "r1", r1, "r1 and r2 cannot both be 0")),
    }
}

fn polygon_corner(sides: u32, radius: f64, flats: bool) -> f64 {
    polygon(sides, radius, flats, 0.0)[0].x
}

/// A wedge on a corner at the origin: its top runs `top_x` along x
/// from the x = 0 side (0: a ramp, built as a triangular prism along
/// y).
pub(super) fn wedge(
    shape: &str,
    [x, y, z]: [f64; 3],
    top_x: f64,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let built = if top_x > 0.0 {
        let corners = |w: f64, h: f64| {
            vec![
                Point3::new(0.0, 0.0, h),
                Point3::new(w, 0.0, h),
                Point3::new(w, y, h),
                Point3::new(0.0, y, h),
            ]
        };
        BRep::polygon_frustum(&corners(x, 0.0), &corners(top_x, z), tolerance)
    } else {
        let section = |at: f64| {
            vec![
                Point3::new(0.0, at, 0.0),
                Point3::new(x, at, 0.0),
                Point3::new(0.0, at, z),
            ]
        };
        BRep::polygon_frustum(&section(0.0), &section(y), tolerance)
    };
    built.map_err(|e| failed(shape, "wedge", e))
}
