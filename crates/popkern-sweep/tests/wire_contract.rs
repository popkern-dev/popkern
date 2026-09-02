//! PK-416: the sweep's wire contract is enforced, not assumed. A
//! profile with a `z` extent (a world-space wire on a plane other than
//! the path frame's) is refused with the measured extent instead of
//! being skinned into collapsed, zero-normal faces — which is what a
//! rectangle lifted onto the YZ plane produced (its `z` spans mapped
//! to nothing because only `x` and `y` reach the frame).

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_sweep::{SweepError, SweepOptions, SweepProfile, sweep};

fn rectangle(points: [Point3; 4]) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: vec![points[0], points[1], points[2], points[3], points[0]],
        weights: vec![1.0; 5],
    }
}

fn path_x() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.5, 0.25), Point3::new(0.6, 0.5, 0.25)],
        weights: vec![1.0, 1.0],
    }
}

#[test]
fn a_wire_with_a_z_extent_is_refused_with_the_measured_extent() {
    // The YZ rectangle as a world-space wire: x = 0, spans along y and z.
    let lifted = rectangle([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.5),
        Point3::new(0.0, 0.0, 0.5),
    ]);
    let error = sweep(
        &SweepProfile::Wire(lifted),
        &path_x(),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect_err("a wire with a z extent must be refused");
    match error {
        SweepError::InvalidProfile(text) => {
            assert!(text.contains("z extent is 5.000e-1"), "{text}");
            assert!(text.contains("path frame"), "{text}");
        }
        other => panic!("{other}"),
    }
}

#[test]
fn the_same_outline_in_the_frame_plane_sweeps_to_a_box_with_live_normals() {
    // Expressed in the path frame (x → normal, y → binormal, z = 0).
    let flat = rectangle([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 0.5, 0.0),
        Point3::new(0.0, 0.5, 0.0),
    ]);
    let result = sweep(
        &SweepProfile::Wire(flat),
        &path_x(),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("served");
    assert_eq!(result.brep.faces.len(), 6);
    assert!(result.warnings.is_empty(), "{:?}", result.warnings);
    for face in result.brep.faces.values() {
        let surface = face.surface.as_ref().expect("surface");
        assert!(
            matches!(
                surface.geometry,
                popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
            ),
            "every face of a swept rectangle is a plane"
        );
    }
    result.brep.validate().expect("valid");
}
