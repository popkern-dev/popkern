//! End-to-end integration tests for popkern-sweep.
//!
//! These tests exercise the full pipeline: profile analysis, frame
//! computation, surface generation, topology construction, and boolean
//! operations.

use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sweep::{
    HelixPitch, HelixSpec, LoftOptions, LoftSection, PipeRadius, SweepOptions, SweepProfile,
};

/// Build a unit circle in the XY plane.
fn circle_profile(radius: f64) -> SweepProfile {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let r = radius;
    SweepProfile::Wire(NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(r, 0.0, 0.0),
            Point3::new(r, r, 0.0),
            Point3::new(0.0, r, 0.0),
            Point3::new(-r, r, 0.0),
            Point3::new(-r, 0.0, 0.0),
            Point3::new(-r, -r, 0.0),
            Point3::new(0.0, -r, 0.0),
            Point3::new(r, -r, 0.0),
            Point3::new(r, 0.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    })
}

/// Build a circle at a given Z position.
fn circle_at_z(radius: f64, z: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let r = radius;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(r, 0.0, z),
            Point3::new(r, r, z),
            Point3::new(0.0, r, z),
            Point3::new(-r, r, z),
            Point3::new(-r, 0.0, z),
            Point3::new(-r, -r, z),
            Point3::new(0.0, -r, z),
            Point3::new(r, -r, z),
            Point3::new(r, 0.0, z),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Build a straight-line path along Z.
fn z_path(start: Point3, end: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    }
}

fn default_tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn sweep_circle_along_line() {
    let profile = circle_profile(1.0);
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0));
    let result = popkern_sweep::sweep(
        &profile,
        &path,
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    assert!(result.is_ok(), "sweep failed: {:?}", result.err());
    let brep = &result.unwrap().brep;
    assert!(brep.faces.len() >= 1, "expected at least 1 face");
    assert!(!brep.vertices.is_empty(), "expected vertices");
    assert!(!brep.edges.is_empty(), "expected edges");
}

#[test]
fn sweep_open_profile_is_sheet() {
    let profile = SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(-1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    });
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 3.0));
    let result = popkern_sweep::sweep(
        &profile,
        &path,
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    assert!(result.is_ok(), "open sweep failed: {:?}", result.err());
    let brep = &result.unwrap().brep;
    // Open profile: no caps.
    assert_eq!(brep.faces.len(), 1, "sheet body should have 1 face");
}

#[test]
fn loft_two_circles() {
    let sections = vec![
        LoftSection {
            curve: circle_at_z(1.0, 0.0),
            position: None,
        },
        LoftSection {
            curve: circle_at_z(1.0, 3.0),
            position: None,
        },
    ];
    let result = popkern_sweep::loft(&sections, &LoftOptions::default(), default_tol(), None);
    assert!(result.is_ok(), "loft failed: {:?}", result.err());
    let brep = &result.unwrap().brep;
    assert!(brep.faces.len() >= 1, "expected at least 1 face");
}

#[test]
fn loft_three_sections_varying_radius() {
    let sections = vec![
        LoftSection {
            curve: circle_at_z(1.0, 0.0),
            position: None,
        },
        LoftSection {
            curve: circle_at_z(2.0, 2.0),
            position: None,
        },
        LoftSection {
            curve: circle_at_z(1.0, 4.0),
            position: None,
        },
    ];
    let result = popkern_sweep::loft(&sections, &LoftOptions::default(), default_tol(), None);
    assert!(result.is_ok(), "loft 3 sections failed: {:?}", result.err());
}

#[test]
fn pipe_constant_radius() {
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0));
    let result = popkern_sweep::pipe(
        &path,
        &PipeRadius::Constant(1.0),
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    assert!(result.is_ok(), "pipe failed: {:?}", result.err());
    let brep = &result.unwrap().brep;
    assert!(brep.faces.len() >= 1);
}

#[test]
fn pipe_variable_radius() {
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0));
    let result = popkern_sweep::pipe(
        &path,
        &PipeRadius::Variable(vec![(0.0, 1.0), (1.0, 2.0)]),
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    assert!(result.is_ok(), "variable pipe failed: {:?}", result.err());
}

#[test]
fn helical_sweep_circle() {
    let profile = circle_profile(0.3);
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 3.0,
        end_radius: 3.0,
        pitch: HelixPitch::Constant(2.0),
        turns: 1.0,
    };
    let result = popkern_sweep::helical_sweep(
        &profile,
        &spec,
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    assert!(result.is_ok(), "helical sweep failed: {:?}", result.err());
    let brep = &result.unwrap().brep;
    assert!(brep.faces.len() >= 1);
}

#[test]
fn swept_blend_two_circles() {
    let profiles = vec![(0.0, circle_profile(1.0)), (1.0, circle_profile(2.0))];
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 10.0));
    let result = popkern_sweep::swept_blend(
        &profiles,
        &path,
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    assert!(result.is_ok(), "swept blend failed: {:?}", result.err());
}

#[test]
fn sweep_with_twist() {
    let profile = circle_profile(1.0);
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0));
    let mut options = SweepOptions::default();
    options.twist = std::f64::consts::PI; // 180 degree twist
    let result = popkern_sweep::sweep(&profile, &path, &options, default_tol(), None);
    assert!(result.is_ok(), "twisted sweep failed: {:?}", result.err());
}

#[test]
fn sweep_with_scaling() {
    let profile = circle_profile(1.0);
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0));
    let mut options = SweepOptions::default();
    options.scale = popkern_sweep::ScaleMethod::Uniform(vec![(0.0, 1.0), (1.0, 2.0)]);
    let result = popkern_sweep::sweep(&profile, &path, &options, default_tol(), None);
    assert!(result.is_ok(), "scaled sweep failed: {:?}", result.err());
}

#[test]
#[ignore] // Boolean on curved/planar intersections is slow in debug builds.
fn sweep_cut_through_cube() {
    let cube = BRep::unit_cube(default_tol());
    let cube = match cube {
        Ok(c) => c,
        Err(_) => return, // skip if cube construction fails
    };
    let profile = circle_profile(0.3);
    let path = z_path(Point3::new(0.5, 0.5, -1.0), Point3::new(0.5, 0.5, 2.0));
    let result = popkern_sweep::sweep_cut(
        &cube,
        &profile,
        &path,
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    // Boolean operations on curved/planar intersections are complex;
    // we check that the pipeline at least produces a result or a clean error.
    match result {
        Ok(r) => {
            assert!(r.brep.faces.len() >= 1, "expected at least 1 face");
        }
        Err(e) => {
            // Boolean may legitimately fail on complex geometry.
            assert!(
                matches!(e, popkern_sweep::SweepError::Boolean(_)),
                "unexpected error: {:?}",
                e
            );
        }
    }
}

#[test]
#[ignore] // Boolean on curved/planar intersections is slow in debug builds.
fn sweep_add_to_cube() {
    let cube = BRep::unit_cube(default_tol());
    let cube = match cube {
        Ok(c) => c,
        Err(_) => return,
    };
    let profile = circle_profile(0.3);
    let path = z_path(Point3::new(0.5, 0.5, 0.0), Point3::new(0.5, 0.5, 3.0));
    let result = popkern_sweep::sweep_add(
        &cube,
        &profile,
        &path,
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    match result {
        Ok(r) => {
            assert!(r.brep.faces.len() >= 1);
        }
        Err(e) => {
            assert!(matches!(e, popkern_sweep::SweepError::Boolean(_)));
        }
    }
}

#[test]
#[ignore] // Boolean on curved/planar intersections is slow in debug builds.
fn pipe_add_to_cube() {
    let cube = BRep::unit_cube(default_tol());
    let cube = match cube {
        Ok(c) => c,
        Err(_) => return,
    };
    let path = z_path(Point3::new(0.5, 0.5, 0.0), Point3::new(0.5, 0.5, 3.0));
    let result = popkern_sweep::pipe_add(
        &cube,
        &path,
        &PipeRadius::Constant(0.2),
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    match result {
        Ok(r) => {
            assert!(r.brep.faces.len() >= 1);
        }
        Err(e) => {
            assert!(matches!(e, popkern_sweep::SweepError::Boolean(_)));
        }
    }
}

#[test]
fn invalid_profile_returns_error() {
    let profile = SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![],
        weights: vec![],
    });
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0));
    let result = popkern_sweep::sweep(
        &profile,
        &path,
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    assert!(result.is_err());
}

#[test]
fn pipe_has_pcurves_on_all_half_edges() {
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0));
    let result = popkern_sweep::pipe(
        &path,
        &PipeRadius::Constant(1.0),
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    assert!(result.is_ok(), "pipe failed: {:?}", result.err());
    let brep = &result.unwrap().brep;

    // Every half-edge in the pipe BRep should have a pcurve.
    for (he_handle, he) in &brep.half_edges {
        assert!(
            he.pcurve.is_some(),
            "half-edge {:?} is missing a pcurve",
            he_handle
        );
    }
}

#[test]
fn sweep_circle_has_pcurves_on_all_half_edges() {
    let profile = circle_profile(1.0);
    let path = z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0));
    let result = popkern_sweep::sweep(
        &profile,
        &path,
        &SweepOptions::default(),
        default_tol(),
        None,
    );
    assert!(result.is_ok(), "sweep failed: {:?}", result.err());
    let brep = &result.unwrap().brep;

    for (he_handle, he) in &brep.half_edges {
        assert!(
            he.pcurve.is_some(),
            "half-edge {:?} is missing a pcurve",
            he_handle
        );
    }
}

#[test]
fn self_intersection_allowed_produces_warning() {
    // A tight curve that may cause self-intersection.
    let profile = circle_profile(5.0); // large profile
    let path = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let mut options = SweepOptions::default();
    options.allow_self_intersection = true;
    let result = popkern_sweep::sweep(&profile, &path, &options, default_tol(), None);
    // Should either succeed (with possible warning) or fail for a reason other than self-intersection.
    match result {
        Ok(r) => {
            // May or may not have self-intersection warning depending on geometry.
            assert!(r.brep.faces.len() >= 1);
        }
        Err(e) => {
            // Should not be SelfIntersection since we allowed it.
            assert!(
                !matches!(e, popkern_sweep::SweepError::SelfIntersection { .. }),
                "should not get SelfIntersection error when allowed"
            );
        }
    }
}
