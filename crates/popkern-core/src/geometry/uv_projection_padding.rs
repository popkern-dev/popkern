//! Scale-aware minimum padding for vertex-bounded parameter windows
//! (PK-462).
//!
//! `bounded_parameter_range` pads a face's projected uv box by half its
//! extent but at least 1.0, and the tight variant by 5% but at least
//! 0.05 — absolute uv lengths (along a cylinder's axis, across a
//! plane). On a millimetre-unit part, whose faces are a few thousandths
//! long, 1.0 is hundreds of times the part: intersection marches
//! traced the extended surfaces across that whole window, leaving the
//! real face and re-projecting at every step (an M6 bolt built 21x
//! slower than the same part at millimetre magnitudes). The floor now
//! scales with the face's own size below 0.5 and is unchanged above,
//! so every model at ordinary scale keeps its windows bit for bit.

use crate::geometry::Point3;
use crate::geometry::surface::SurfaceGeometry;

/// `floor`, scaled by the face's size when that is below 0.5: the
/// larger of the vertices' 3D extent and the surface's own radius (a
/// cap with one vertex still has its circle). Degenerate sizes keep
/// the legacy floor.
pub(crate) fn scaled_padding_floor(
    floor: f64,
    vertices: &[Point3],
    surface: &SurfaceGeometry,
) -> f64 {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in vertices {
        for (i, c) in [p.x, p.y, p.z].into_iter().enumerate() {
            lo[i] = lo[i].min(c);
            hi[i] = hi[i].max(c);
        }
    }
    let extent = if vertices.is_empty() {
        0.0
    } else {
        ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt()
    };
    let radius = match surface {
        SurfaceGeometry::Cylinder { radius, .. } | SurfaceGeometry::Sphere { radius, .. } => {
            *radius
        }
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } => major_radius + minor_radius,
        _ => 0.0,
    };
    let size = extent.max(radius.abs());
    if !(size.is_finite() && size > 0.0) {
        return floor;
    }
    floor * (2.0 * size).min(1.0)
}

#[cfg(test)]
mod tests {
    use super::scaled_padding_floor;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::{Point3, Vec3};

    fn cylinder(radius: f64) -> SurfaceGeometry {
        SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius,
        }
    }

    #[test]
    fn ordinary_scale_faces_keep_the_legacy_floor() {
        let v = [Point3::new(3.0, 0.0, 0.0), Point3::new(3.0, 0.0, 10.0)];
        assert_eq!(scaled_padding_floor(1.0, &v, &cylinder(3.0)), 1.0);
        assert_eq!(scaled_padding_floor(0.05, &v, &cylinder(3.0)), 0.05);
        // A unit cube's faces are at the boundary: unchanged.
        let q = [Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        assert_eq!(scaled_padding_floor(1.0, &q, &cylinder(0.0)), 1.0);
    }

    #[test]
    fn millimetre_scale_faces_get_a_proportional_floor() {
        let s = 1e-3;
        let v = [
            Point3::new(3.0 * s, 0.0, 0.0),
            Point3::new(3.0 * s, 0.0, 10.0 * s),
        ];
        let floor = scaled_padding_floor(1.0, &v, &cylinder(3.0 * s));
        let extent = (10.0f64 * s).hypot(0.0);
        assert!((floor - 2.0 * extent).abs() < 1e-15, "floor {floor}");
    }

    #[test]
    fn degenerate_faces_keep_the_legacy_floor() {
        assert_eq!(scaled_padding_floor(1.0, &[], &cylinder(0.0)), 1.0);
    }
}
