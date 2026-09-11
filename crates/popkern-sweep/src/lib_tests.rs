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
        start_angle: 0.0,
        left_handed: false,
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
    let count = compute_sample_count(&table, &path);
    assert!(count >= MIN_FRAMES);
    assert!(count <= 1024);
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

/// PK-477: a genuine multi-wire closed sweep (outer + inner square)
/// builds a valid genus-1 solid — per-wire rings, caps with inner
/// loops. Before the per-wire builder the flattened ring
/// cross-stitched the wires (measured: validate() failed with a
/// 1.414 loop gap).
#[test]
fn multi_wire_square_tube_sweeps_valid() {
    let sq = |h: f64| square_wire(h);
    let path = z_path(5.0);
    let profile = SweepProfile::MultiWire(vec![sq(2.0), sq(1.0)]);
    let result = sweep(
        &profile,
        &path,
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("multi-wire sweep");
    let brep = &result.brep;
    brep.validate().expect("valid tube solid");
    // 4 outer walls + 4 inner walls + 2 caps.
    assert_eq!(brep.faces.len(), 10, "tube face census");
    // Genus-1: V - E + 2F - L == 0.
    let chi = brep.vertices.len() as i64 - brep.edges.len() as i64 + 2 * brep.faces.len() as i64
        - brep.loops.len() as i64;
    assert_eq!(chi, 0, "tube is genus-1");
}

/// PK-477: multi-wire winding is normalized per role regardless of
/// authored direction — a CW outer and CCW inner produce the same
/// valid tube.
#[test]
fn multi_wire_winding_is_normalized() {
    let mut outer = square_wire(2.0);
    outer.control_points.reverse();
    outer.weights.reverse();
    let inner = square_wire(1.0);
    let path = z_path(5.0);
    let profile = SweepProfile::MultiWire(vec![outer, inner]);
    let result = sweep(
        &profile,
        &path,
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("multi-wire sweep");
    result.brep.validate().expect("valid tube solid");
}

/// PK-477: honest declines — an open wire in a multi-wire profile,
/// and an "inner" wire lying outside the outer.
#[test]
fn multi_wire_invalid_profiles_decline() {
    let path = z_path(5.0);
    let open_wire = {
        let mut w = square_wire(1.0);
        w.control_points.pop();
        w.weights.pop();
        let n = w.control_points.len();
        w.knots = std::iter::once(0.0)
            .chain((0..n).map(|i| i as f64 / (n - 1) as f64))
            .chain(std::iter::once(1.0))
            .collect();
        w
    };
    let r = sweep(
        &SweepProfile::MultiWire(vec![square_wire(2.0), open_wire]),
        &path,
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    );
    assert!(r.is_err(), "open inner wire must decline");
    let mut outside = square_wire(1.0);
    for p in &mut outside.control_points {
        p.x += 10.0;
    }
    let r = sweep(
        &SweepProfile::MultiWire(vec![square_wire(2.0), outside]),
        &path,
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    );
    assert!(r.is_err(), "disjoint inner wire must decline");
}

/// Degree-1 closed square wire of half-width `h` in the XY plane.
fn square_wire(h: f64) -> NurbsCurve {
    let pts = [[-h, -h], [h, -h], [h, h], [-h, h], [-h, -h]];
    let n = pts.len();
    let mut knots = vec![0.0];
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    NurbsCurve {
        degree: 1,
        knots,
        control_points: pts
            .iter()
            .map(|p| popkern_core::geometry::Point3::new(p[0], p[1], 0.0))
            .collect(),
        weights: vec![1.0; n],
    }
}

/// PK-477: the multi-wire path composes with a helical path (frame
/// machinery smoke — the per-wire rings are path-agnostic).
#[test]
fn multi_wire_helical_sweep_smoke() {
    use crate::{FrameMethod, HelixPitch, HelixSpec, helical_sweep};
    let profile = SweepProfile::MultiWire(vec![square_wire(0.5), square_wire(0.2)]);
    let helix = HelixSpec {
        axis_origin: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        axis_dir: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        start_radius: 5.0,
        end_radius: 5.0,
        pitch: HelixPitch::Constant(3.0),
        turns: 0.75,
        start_angle: 0.0,
        left_handed: false,
    };
    let options = SweepOptions {
        frame: FrameMethod::Frenet,
        ..SweepOptions::default()
    };
    let result = helical_sweep(&profile, &helix, &options, Tolerance::default(), None)
        .expect("multi-wire helical sweep");
    result.brep.validate().expect("valid helical tube");
}
