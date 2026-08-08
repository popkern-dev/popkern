//! PK-245: winding-flag honesty audit at the sweep serving boundary.
//!
//! Every consumer of the stored orientation flag (edge convexity
//! classification, mesher triangle winding, STEP `same_sense`,
//! measured volume sign) rests on the interior-left invariant the
//! core audit measures: outer loop CCW in UV ⇔ flag true. PK-244
//! closed the invariant at the fillet/chamfer/boolean tails; this
//! file measures the SWEEP operator family, one fixture per entry
//! point. A failure lists the dishonest faces with their surface
//! class — that IS the defect record for the author behind them.

use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sweep::{
    HelixPitch, HelixSpec, LoftOptions, LoftSection, PipeRadius, SweepOptions, SweepProfile,
};

/// Format every dishonest face with its surface class for the
/// assertion message; empty = the operator's authors are honest.
fn dishonest_report(brep: &BRep) -> Vec<String> {
    brep.audit_winding_honesty()
        .into_iter()
        .filter(|v| v.honest == Some(false))
        .map(|v| {
            let surface = brep
                .faces
                .get(v.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)))
                .unwrap_or_else(|| "<no-surface>".into());
            format!(
                "{:?} [{surface}] area {:.3e} flag {}",
                v.face, v.signed_uv_area, v.orientation
            )
        })
        .collect()
}

fn assert_honest(brep: &BRep, fixture: &str) {
    let report = dishonest_report(brep);
    assert!(
        report.is_empty(),
        "{fixture}: {} dishonest face(s):\n  {}",
        report.len(),
        report.join("\n  ")
    );
}

fn circle_profile(radius: f64) -> SweepProfile {
    SweepProfile::Wire(circle_at_z(radius, 0.0))
}

fn circle_at_z(r: f64, z: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
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

fn z_path(start: Point3, end: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    }
}

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Closed profile swept along a line: side surface + two caps.
#[test]
fn sweep_solid_is_winding_honest() {
    let result = popkern_sweep::sweep(
        &circle_profile(1.0),
        &z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 5.0)),
        &SweepOptions::default(),
        tol(),
        None,
    )
    .expect("sweep");
    assert_honest(&result.brep, "sweep_solid");
}

/// Open profile swept along a line: single-face sheet body.
#[test]
fn sweep_open_sheet_is_winding_honest() {
    let profile = SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(-1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    });
    let result = popkern_sweep::sweep(
        &profile,
        &z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 3.0)),
        &SweepOptions::default(),
        tol(),
        None,
    )
    .expect("open sweep");
    assert_honest(&result.brep, "sweep_open_sheet");
}

/// Loft through two circular sections.
#[test]
fn loft_solid_is_winding_honest() {
    let sections = vec![
        LoftSection {
            curve: circle_at_z(1.0, 0.0),
            position: None,
        },
        LoftSection {
            curve: circle_at_z(0.6, 2.0),
            position: None,
        },
    ];
    let result =
        popkern_sweep::loft(&sections, &LoftOptions::default(), tol(), None).expect("loft");
    assert_honest(&result.brep, "loft_solid");
}

/// Constant-radius pipe along a straight path.
#[test]
fn pipe_solid_is_winding_honest() {
    let result = popkern_sweep::pipe(
        &z_path(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 4.0)),
        &PipeRadius::Constant(0.5),
        &SweepOptions::default(),
        tol(),
        None,
    )
    .expect("pipe");
    assert_honest(&result.brep, "pipe_solid");
}

/// Helical sweep: curved path exercises frame-driven surface authors.
#[test]
fn helical_sweep_is_winding_honest() {
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 2.0,
        end_radius: 2.0,
        pitch: HelixPitch::Constant(1.0),
        turns: 2.0,
    };
    let result = popkern_sweep::helical_sweep(
        &circle_profile(0.3),
        &spec,
        &SweepOptions::default(),
        tol(),
        None,
    )
    .expect("helical sweep");
    assert_honest(&result.brep, "helical_sweep");
}
