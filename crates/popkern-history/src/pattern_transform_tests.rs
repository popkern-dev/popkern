// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//! Tests for per-instance transform computation.
#[cfg(test)]
use super::*;

use crate::feature::{AxisDesc, Feature, Vec3Desc};
use crate::parameter::ParameterTable;
use crate::pattern::{MirrorPlaneDesc, PatternDirection};
use crate::unit::Unit;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Build a minimal [`ParameterTable`] and return two parameter IDs:
/// `(count_id, dist_or_angle_id)`.
fn make_params(
    count: f64,
    dist: f64,
) -> (
    ParameterTable,
    crate::parameter::ParameterId,
    crate::parameter::ParameterId,
) {
    let mut table = ParameterTable::new();
    let count_id = table.add_global("count", count, Unit::Ratio).unwrap();
    let dist_id = table.add_global("dist", dist, Unit::Ratio).unwrap();
    (table, count_id, dist_id)
}

// ---------------------------------------------------------------------------
// Linear pattern
// ---------------------------------------------------------------------------

#[test]
fn linear_4_instances_along_x() {
    let brep = BRep::new(Tolerance::default());
    let (table, count_id, dist_id) = make_params(4.0, 30.0);

    let feature = Feature::LinearPattern {
        seed_features: vec![],
        primary: PatternDirection {
            direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
            count: count_id,
            total_distance: dist_id,
        },
        secondary: None,
        geometry_pattern: false,
        instance_overrides: Default::default(),
    };

    let transforms = compute_pattern_transforms(&feature, &brep, &table).unwrap();
    assert_eq!(transforms.len(), 4);

    let expected = [0.0, 10.0, 20.0, 30.0];
    for (i, &x) in expected.iter().enumerate() {
        let t = &transforms[i];
        // translation along X: matrix[0][3]
        assert!(
            (t.matrix[0][3] - x).abs() < 1e-10,
            "instance {i}: expected tx={x}, got {}",
            t.matrix[0][3]
        );
        assert!(t.matrix[1][3].abs() < 1e-10, "instance {i}: ty should be 0");
        assert!(t.matrix[2][3].abs() < 1e-10, "instance {i}: tz should be 0");
    }
}

#[test]
fn linear_1_instance_is_identity() {
    let brep = BRep::new(Tolerance::default());
    let (table, count_id, dist_id) = make_params(1.0, 10.0);

    let feature = Feature::LinearPattern {
        seed_features: vec![],
        primary: PatternDirection {
            direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
            count: count_id,
            total_distance: dist_id,
        },
        secondary: None,
        geometry_pattern: false,
        instance_overrides: Default::default(),
    };

    let transforms = compute_pattern_transforms(&feature, &brep, &table).unwrap();
    assert_eq!(transforms.len(), 1);
    assert!(transforms[0].approx_eq(&Transform3::identity(), 1e-10));
}

#[test]
fn linear_2d_grid_3x2() {
    let brep = BRep::new(Tolerance::default());
    let mut table = ParameterTable::new();
    let count1_id = table.add_global("c1", 3.0, Unit::Ratio).unwrap();
    let dist1_id = table.add_global("d1", 20.0, Unit::Ratio).unwrap();
    let count2_id = table.add_global("c2", 2.0, Unit::Ratio).unwrap();
    let dist2_id = table.add_global("d2", 10.0, Unit::Ratio).unwrap();

    let feature = Feature::LinearPattern {
        seed_features: vec![],
        primary: PatternDirection {
            direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
            count: count1_id,
            total_distance: dist1_id,
        },
        secondary: Some(PatternDirection {
            direction: Vec3Desc::Literal([0.0, 1.0, 0.0]),
            count: count2_id,
            total_distance: dist2_id,
        }),
        geometry_pattern: false,
        instance_overrides: Default::default(),
    };

    let transforms = compute_pattern_transforms(&feature, &brep, &table).unwrap();
    assert_eq!(transforms.len(), 6, "3×2 grid should have 6 instances");

    // Row-major ordering: index = i * count2 + j.
    // i=0, j=0 → (0,0);  i=0, j=1 → (0,10)
    // i=1, j=0 → (10,0); i=1, j=1 → (10,10)
    // i=2, j=0 → (20,0); i=2, j=1 → (20,10)
    let expected = [
        (0.0, 0.0),
        (0.0, 10.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (20.0, 0.0),
        (20.0, 10.0),
    ];
    for (idx, (ex, ey)) in expected.iter().enumerate() {
        let t = &transforms[idx];
        assert!(
            (t.matrix[0][3] - ex).abs() < 1e-10,
            "instance {idx}: expected tx={ex}, got {}",
            t.matrix[0][3]
        );
        assert!(
            (t.matrix[1][3] - ey).abs() < 1e-10,
            "instance {idx}: expected ty={ey}, got {}",
            t.matrix[1][3]
        );
    }
}

#[test]
fn linear_count_zero_returns_error() {
    let brep = BRep::new(Tolerance::default());
    let (table, count_id, dist_id) = make_params(0.0, 10.0);

    let feature = Feature::LinearPattern {
        seed_features: vec![],
        primary: PatternDirection {
            direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
            count: count_id,
            total_distance: dist_id,
        },
        secondary: None,
        geometry_pattern: false,
        instance_overrides: Default::default(),
    };

    let result = compute_pattern_transforms(&feature, &brep, &table);
    assert!(result.is_err(), "count=0 should return an error");
}

// ---------------------------------------------------------------------------
// Circular pattern
// ---------------------------------------------------------------------------

#[test]
fn circular_4_instances_full_circle() {
    // With count=4 and total_angle = 3π/2, step = total_angle/(count-1) = π/2.
    // Instances land at 0°, 90°, 180°, 270°.
    let brep = BRep::new(Tolerance::default());
    let mut table = ParameterTable::new();
    let count_id = table.add_global("count", 4.0, Unit::Ratio).unwrap();
    // total_angle = (count-1) * π/2 = 3π/2 gives steps of exactly π/2.
    let total_angle = 3.0 * std::f64::consts::FRAC_PI_2;
    let angle_id = table.add_global("angle", total_angle, Unit::Ratio).unwrap();

    let feature = Feature::CircularPattern {
        seed_features: vec![],
        axis: AxisDesc::Literal {
            origin: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
        },
        count: count_id,
        total_angle: angle_id,
        geometry_pattern: false,
        instance_overrides: Default::default(),
    };

    let transforms = compute_pattern_transforms(&feature, &brep, &table).unwrap();
    assert_eq!(transforms.len(), 4);

    // Instance 0 is identity.
    assert!(
        transforms[0].approx_eq(&Transform3::identity(), 1e-10),
        "instance 0 should be identity"
    );

    // Steps of π/2: instances at 0, π/2, π, 3π/2.
    let expected_angles = [
        0.0,
        std::f64::consts::FRAC_PI_2,
        std::f64::consts::PI,
        3.0 * std::f64::consts::FRAC_PI_2,
    ];
    for (i, &expected_angle) in expected_angles.iter().enumerate() {
        let ref_transform = Transform3::rotation_around_axis(
            popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
            expected_angle,
        )
        .unwrap();
        assert!(
            transforms[i].approx_eq(&ref_transform, 1e-10),
            "instance {i}: rotation mismatch for angle {expected_angle}"
        );
    }
}

#[test]
fn circular_1_instance_is_identity() {
    let brep = BRep::new(Tolerance::default());
    let mut table = ParameterTable::new();
    let count_id = table.add_global("count", 1.0, Unit::Ratio).unwrap();
    let angle_id = table
        .add_global("angle", std::f64::consts::TAU, Unit::Ratio)
        .unwrap();

    let feature = Feature::CircularPattern {
        seed_features: vec![],
        axis: AxisDesc::Literal {
            origin: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
        },
        count: count_id,
        total_angle: angle_id,
        geometry_pattern: false,
        instance_overrides: Default::default(),
    };

    let transforms = compute_pattern_transforms(&feature, &brep, &table).unwrap();
    assert_eq!(transforms.len(), 1);
    assert!(
        transforms[0].approx_eq(&Transform3::identity(), 1e-10),
        "single instance should be identity"
    );
}

#[test]
fn circular_count_zero_returns_error() {
    let brep = BRep::new(Tolerance::default());
    let mut table = ParameterTable::new();
    let count_id = table.add_global("count", 0.0, Unit::Ratio).unwrap();
    let angle_id = table.add_global("angle", 1.0, Unit::Ratio).unwrap();

    let feature = Feature::CircularPattern {
        seed_features: vec![],
        axis: AxisDesc::Literal {
            origin: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
        },
        count: count_id,
        total_angle: angle_id,
        geometry_pattern: false,
        instance_overrides: Default::default(),
    };

    let result = compute_pattern_transforms(&feature, &brep, &table);
    assert!(result.is_err(), "count=0 should return an error");
}

// ---------------------------------------------------------------------------
// Mirror pattern
// ---------------------------------------------------------------------------

#[test]
fn mirror_yz_gives_identity_and_reflection() {
    let brep = BRep::new(Tolerance::default());

    let feature = Feature::Mirror {
        seed_features: vec![],
        plane: MirrorPlaneDesc::YZ,
        merge_result: true,
    };

    let transforms = compute_pattern_transforms(&feature, &brep, &ParameterTable::new()).unwrap();
    assert_eq!(transforms.len(), 2);

    // First is identity.
    assert!(
        transforms[0].approx_eq(&Transform3::identity(), 1e-10),
        "first mirror transform should be identity"
    );

    // Second is reflection across the YZ plane (normal = X axis).
    let expected_reflection = Transform3::reflection(
        popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0),
    )
    .unwrap();
    assert!(
        transforms[1].approx_eq(&expected_reflection, 1e-10),
        "second mirror transform should be YZ-plane reflection"
    );
}

#[test]
fn mirror_xy_gives_correct_reflection() {
    let brep = BRep::new(Tolerance::default());

    let feature = Feature::Mirror {
        seed_features: vec![],
        plane: MirrorPlaneDesc::XY,
        merge_result: true,
    };

    let transforms = compute_pattern_transforms(&feature, &brep, &ParameterTable::new()).unwrap();
    assert_eq!(transforms.len(), 2);

    let expected_reflection = Transform3::reflection(
        popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
    )
    .unwrap();
    assert!(
        transforms[1].approx_eq(&expected_reflection, 1e-10),
        "second mirror transform should be XY-plane reflection"
    );
}

#[test]
fn mirror_custom_plane() {
    let brep = BRep::new(Tolerance::default());

    let feature = Feature::Mirror {
        seed_features: vec![],
        plane: MirrorPlaneDesc::Custom {
            origin: [1.0, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
        },
        merge_result: false,
    };

    let transforms = compute_pattern_transforms(&feature, &brep, &ParameterTable::new()).unwrap();
    assert_eq!(transforms.len(), 2);
    assert!(
        transforms[0].approx_eq(&Transform3::identity(), 1e-10),
        "first transform should be identity"
    );

    let expected_reflection = Transform3::reflection(
        popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
        popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0),
    )
    .unwrap();
    assert!(
        transforms[1].approx_eq(&expected_reflection, 1e-10),
        "reflection mismatch for custom plane"
    );
}

#[test]
fn non_pattern_feature_returns_error() {
    let brep = BRep::new(Tolerance::default());
    let feature = Feature::Sketch {
        plane: crate::feature::SketchPlaneDesc::XY,
        sketch_data: String::new(),
    };
    let result = compute_pattern_transforms(&feature, &brep, &ParameterTable::new());
    assert!(
        result.is_err(),
        "non-pattern feature should return an error"
    );
}
