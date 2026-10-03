//! Axis completion for placements whose reference direction is omitted.

use popkern_core::geometry::Vec3;

/// The x axis of an `AXIS2_PLACEMENT_3D` whose `ref_direction` is
/// omitted (ISO 10303-42 `first_proj_axis`): +X projected orthogonal to
/// the axis, or +Y when the axis is along X. A constant +X default made
/// every X-normal plane's frame degenerate (measured: `clevis21.stp`'s
/// X-normal planes mapped their in-plane edges 29.5 units off the
/// surface and their faces meshed open).
pub(super) fn default_ref_direction(z_axis: Vec3) -> Vec3 {
    let x = Vec3::new(1.0, 0.0, 0.0);
    let Some(z) = z_axis.normalized() else {
        return x;
    };
    let arbitrary = if z.cross(&x).length() < 1e-9 {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        x
    };
    (arbitrary - z * arbitrary.dot(&z))
        .normalized()
        .unwrap_or(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_is_orthogonal_to_every_axis() {
        for z in [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 1.0, 1.0),
        ] {
            let x = default_ref_direction(z);
            assert!((x.length() - 1.0).abs() < 1e-12);
            assert!(
                x.dot(&z.normalized().expect("unit")).abs() < 1e-12,
                "{z:?} -> {x:?}"
            );
        }
        // The common case is the spec's +X.
        let x = default_ref_direction(Vec3::new(0.0, 0.0, 1.0));
        assert!((x.x - 1.0).abs() < 1e-15);
    }
}
