// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::*;
use popkern_core::geometry::Point3;
use popkern_core::history::TrackingContext;
use popkern_core::naming::OperationId;

/// Build a unit circle profile.
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

/// Build a straight-line path along the Z axis.
fn z_path(length: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, length)],
        weights: vec![1.0, 1.0],
    }
}

/// Build a line segment profile (open).
fn line_profile() -> SweepProfile {
    SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(-1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    })
}

#[test]
fn sweep_circle_along_z() {
    let profile = circle_profile(1.0);
    let path = z_path(5.0);
    let result = sweep(
        &profile,
        &path,
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "sweep failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    assert!(!brep.faces.is_empty(), "expected at least 1 face");
}

#[test]
fn sweep_open_profile_produces_sheet() {
    let profile = line_profile();
    let path = z_path(3.0);
    let result = sweep(
        &profile,
        &path,
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "sweep failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    // Open profile: no caps, just side faces.
    assert!(!brep.faces.is_empty());
}

#[test]
fn loft_two_circles() {
    let sections = vec![
        LoftSection {
            curve: circle_profile(1.0).as_wire(),
            position: None,
        },
        LoftSection {
            curve: circle_profile(1.0).as_wire_at_z(3.0),
            position: None,
        },
    ];
    let result = loft(
        &sections,
        &LoftOptions::default(),
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "loft failed: {:?}", result.err());
}

#[test]
fn pipe_constant_radius() {
    let path = z_path(5.0);
    let result = pipe(
        &path,
        &PipeRadius::Constant(1.0),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "pipe failed: {:?}", result.err());
}

#[test]
fn helical_sweep_circle() {
    let profile = circle_profile(0.5);
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        start_radius: 3.0,
        end_radius: 3.0,
        pitch: HelixPitch::Constant(2.0),
        turns: 1.0,
    };
    let result = helical_sweep(
        &profile,
        &spec,
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "helical sweep failed: {:?}", result.err());
}

#[test]
fn compute_sample_count_reasonable() {
    let path = z_path(5.0);
    let table = arc_length::ArcLengthTable::build(&path, 100);
    let count = compute_sample_count(&table);
    assert!(count >= MIN_FRAMES);
    assert!(count <= 200);
}

#[test]
fn check_path_closed_open() {
    let path = z_path(5.0);
    let tol = Tolerance::default();
    assert!(!check_path_closed(&path, &tol));
}

#[test]
fn sweep_tracking_assigns_names() {
    let profile = circle_profile(1.0);
    let path = z_path(5.0);
    let mut ctx = TrackingContext::new(OperationId(42));
    let result = sweep(
        &profile,
        &path,
        &SweepOptions::default(),
        Tolerance::default(),
        Some(&mut ctx),
    );
    assert!(result.is_ok(), "sweep failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests

    // All faces should have names assigned.
    for (_, face) in &brep.faces {
        assert!(face.name.is_some(), "face should have a name");
    }
    // All edges should have names assigned.
    for (_, edge) in &brep.edges {
        assert!(edge.name.is_some(), "edge should have a name");
    }
    // All vertices should have names assigned.
    for (_, vertex) in &brep.vertices {
        assert!(vertex.name.is_some(), "vertex should have a name");
    }

    // The tracking context should have recorded generated entities.
    let history = ctx.finish();
    assert!(
        !history.generated.is_empty(),
        "tracking should record generated entities"
    );
    // Total generated = faces + edges + vertices.
    let expected = brep.faces.len() + brep.edges.len() + brep.vertices.len();
    assert_eq!(history.generated.len(), expected);
}

// Helper trait for test convenience.
impl SweepProfile {
    /// Extract the inner wire curve (panics if not a Wire variant).
    fn as_wire(&self) -> NurbsCurve {
        match self {
            SweepProfile::Wire(c) => c.clone(),
            _ => panic!("expected Wire"),
        }
    }

    /// Extract the wire and translate to a given Z position.
    fn as_wire_at_z(&self, z: f64) -> NurbsCurve {
        match self {
            SweepProfile::Wire(c) => {
                let control_points = c
                    .control_points
                    .iter()
                    .map(|p| popkern_core::geometry::Point3::new(p.x, p.y, z))
                    .collect();
                NurbsCurve {
                    degree: c.degree,
                    knots: c.knots.clone(),
                    control_points,
                    weights: c.weights.clone(),
                }
            }
            _ => panic!("expected Wire"),
        }
    }
}
