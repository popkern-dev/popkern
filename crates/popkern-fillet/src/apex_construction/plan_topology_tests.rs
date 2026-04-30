//! Unit tests for [`decide_cap_kind`].
//!
//! Tests cap-kind dispatch (Sphere vs Gregory) on synthetic
//! [`PerFilletTermination`] inputs, covering ADR §4.1's tolerance-
//! aware criteria.

use super::cap_kind::CapKind;
use super::plan_terminations::PerFilletTermination;
use super::plan_topology::decide_cap_kind;
use super::telemetry::PlanTelemetry;
use crate::types::{FilletRadius, FilletSpec};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{EdgeHandle, FaceHandle};

fn dummy_curve() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

fn term(
    idx: usize,
    endpoint_a: Point3,
    endpoint_b: Point3,
    cap_centre: Point3,
    cap_radius: f64,
) -> PerFilletTermination {
    PerFilletTermination {
        fillet_idx: idx,
        endpoint_a,
        endpoint_b,
        cap_centre,
        cap_radius,
        curve_3d: dummy_curve(),
        face_a_normal: Vec3::new(1.0, 0.0, 0.0),
        face_b_normal: Vec3::new(0.0, 1.0, 0.0),
        face_a: FaceHandle::default(),
        face_b: FaceHandle::default(),
    }
}

fn dummy_specs(n: usize, r: f64) -> Vec<FilletSpec> {
    (0..n)
        .map(|_| FilletSpec {
            edge: EdgeHandle::default(),
            radius: FilletRadius::Constant(r),
        })
        .collect()
}

#[test]
fn dispatches_sphere_for_uniform_radii_and_collapsed_bridges() {
    // Three fillets at uniform radius 0.1, bridges collapse
    // (T_b[k] == T_a[k+1] within tol).
    let p0 = Point3::new(0.1, 0.0, 0.1);
    let p1 = Point3::new(0.1, 0.1, 0.0);
    let p2 = Point3::new(0.0, 0.1, 0.1);
    let centre = Point3::new(0.0, 0.0, 0.0);
    let terms = vec![
        term(0, p2, p0, centre, 0.1),
        term(1, p0, p1, centre, 0.1),
        term(2, p1, p2, centre, 0.1),
    ];
    let specs = dummy_specs(3, 0.1);
    let mut telemetry = PlanTelemetry::new();
    let kind = decide_cap_kind(&terms, &specs, &Tolerance::default(), &mut telemetry)
        .expect("decide_cap_kind should succeed");
    assert!(
        matches!(kind, CapKind::Sphere { .. }),
        "expected Sphere, got {kind:?}"
    );
    assert!(telemetry.cap_kind_reason.contains("Sphere"));
}

#[test]
fn dispatches_gregory_for_variable_radii() {
    // Three fillets, radii 0.08 / 0.10 / 0.12 (not uniform).
    let p0 = Point3::new(0.1, 0.0, 0.1);
    let p1 = Point3::new(0.1, 0.1, 0.0);
    let p2 = Point3::new(0.0, 0.1, 0.1);
    let centre = Point3::new(0.0, 0.0, 0.0);
    let terms = vec![
        term(0, p2, p0, centre, 0.08),
        term(1, p0, p1, centre, 0.10),
        term(2, p1, p2, centre, 0.12),
    ];
    let specs = dummy_specs(3, 0.10);
    let mut telemetry = PlanTelemetry::new();
    let kind = decide_cap_kind(&terms, &specs, &Tolerance::default(), &mut telemetry)
        .expect("decide_cap_kind should succeed");
    assert!(
        matches!(kind, CapKind::Gregory { .. }),
        "expected Gregory, got {kind:?}"
    );
    assert!(telemetry.cap_kind_reason.contains("variable radii"));
}

#[test]
fn dispatches_gregory_for_uniform_radii_but_disjoint_bridges() {
    // Uniform radii but T_b[k] != T_a[k+1] (bridges don't collapse).
    let p0a = Point3::new(0.1, 0.0, 0.1);
    let p0b = Point3::new(0.1, 0.05, 0.1);
    let p1a = Point3::new(0.1, 0.1, 0.0);
    let p1b = Point3::new(0.05, 0.1, 0.05);
    let p2a = Point3::new(0.0, 0.1, 0.1);
    let p2b = Point3::new(0.05, 0.05, 0.1);
    let centre = Point3::new(0.0, 0.0, 0.0);
    let terms = vec![
        term(0, p0a, p0b, centre, 0.1),
        term(1, p1a, p1b, centre, 0.1),
        term(2, p2a, p2b, centre, 0.1),
    ];
    let specs = dummy_specs(3, 0.1);
    let mut telemetry = PlanTelemetry::new();
    let kind = decide_cap_kind(&terms, &specs, &Tolerance::default(), &mut telemetry)
        .expect("decide_cap_kind should succeed");
    assert!(
        matches!(kind, CapKind::Gregory { .. }),
        "expected Gregory (bridges open), got {kind:?}"
    );
}

#[test]
fn telemetry_records_diagnostic_values() {
    let p0 = Point3::new(0.1, 0.0, 0.1);
    let p1 = Point3::new(0.1, 0.1, 0.0);
    let p2 = Point3::new(0.0, 0.1, 0.1);
    let centre = Point3::new(0.0, 0.0, 0.0);
    let terms = vec![
        term(0, p2, p0, centre, 0.08),
        term(1, p0, p1, centre, 0.10),
        term(2, p1, p2, centre, 0.12),
    ];
    let specs = dummy_specs(3, 0.10);
    let mut telemetry = PlanTelemetry::new();
    let _ = decide_cap_kind(&terms, &specs, &Tolerance::default(), &mut telemetry).unwrap();
    assert!(telemetry.diagnostic_values.contains_key("radii_min"));
    assert!(telemetry.diagnostic_values.contains_key("radii_max"));
    assert!(telemetry.diagnostic_values.contains_key("radii_spread"));
    assert!(
        telemetry
            .diagnostic_values
            .contains_key("max_bridge_collapse")
    );
    assert!(
        telemetry
            .diagnostic_values
            .contains_key("bridge_collapse_threshold")
    );
}
