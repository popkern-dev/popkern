use super::split_pcurve_at_image;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};

fn cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    }
}

#[test]
fn linear_pcurve_splits_at_the_image() {
    use crate::geometry::surface::SurfaceEval;
    // Vertical iso line on the cylinder at theta = 0: (0, z) chart.
    let pc = Pcurve::linear((0.0, 0.0), (0.0, 2.0), 0.0, 1.0);
    let surface = cylinder();
    let split = surface.point_at(0.0, 0.5);
    let (a, b) = split_pcurve_at_image(&pc, &surface, &split, 1.0e-6).expect("split");
    let a_end = a.point_at(a.t_end);
    let b_start = b.point_at(b.t_start);
    assert!((a_end.y - 0.5).abs() < 1.0e-9, "cut v: {}", a_end.y);
    assert!((b_start.y - 0.5).abs() < 1.0e-9);
    // Branch preserved: both pieces stay at u = 0.
    assert!(a.point_at(a.t_start).x.abs() < 1.0e-12);
    assert!(b.point_at(b.t_end).x.abs() < 1.0e-12);
}

#[test]
fn split_declines_when_the_image_never_reaches_the_point() {
    use crate::geometry::surface::SurfaceEval;
    let pc = Pcurve::linear((0.0, 0.0), (0.0, 2.0), 0.0, 1.0);
    let surface = cylinder();
    // A point on the far side of the cylinder, off the iso image.
    let split = surface.point_at(std::f64::consts::PI, 0.5);
    assert!(split_pcurve_at_image(&pc, &surface, &split, 1.0e-6).is_none());
}

#[test]
fn split_declines_at_an_endpoint_cut() {
    use crate::geometry::surface::SurfaceEval;
    let pc = Pcurve::linear((0.0, 0.0), (0.0, 2.0), 0.0, 1.0);
    let surface = cylinder();
    let split = surface.point_at(0.0, 0.0);
    assert!(split_pcurve_at_image(&pc, &surface, &split, 1.0e-6).is_none());
}

#[test]
fn fraction_split_covers_the_whole_range() {
    let pc = Pcurve::linear((0.2, 0.0), (0.8, 1.0), 0.0, 4.0);
    let (a, b) = pc.split_at_fraction(0.25).expect("split");
    assert!((a.t_start - 0.0).abs() < 1e-12);
    assert!((a.t_end - 1.0).abs() < 1e-12);
    assert!((b.t_start - 1.0).abs() < 1e-12);
    assert!((b.t_end - 4.0).abs() < 1e-12);
    let mid_a = a.point_at(a.t_end);
    let mid_b = b.point_at(b.t_start);
    assert!((mid_a.x - mid_b.x).abs() < 1e-12);
    assert!((mid_a.y - mid_b.y).abs() < 1e-12);
}

/// PK-669: a pcurve whose image WINDS — the helical root edge of an
/// M2.5 thread tool, 13.9 turns at radius 0.98 and pitch 0.45 — has
/// one distance basin per turn. A 64-sample scan spaces samples 78°
/// apart, an adjacent turn's sample (one pitch away) can beat the
/// true turn's, and the old global-best refinement converged a pitch
/// off and declined the split. Every split point along the winding
/// must cut exactly at its own image.
#[test]
fn winding_image_splits_on_the_true_turn() {
    use crate::geometry::surface::SurfaceEval;
    let radius = 0.98;
    let pitch = 0.45;
    let turns = 13.9;
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    };
    let tau = std::f64::consts::TAU;
    let pc = Pcurve::linear((0.0, 0.0), (turns * tau, turns * pitch), 0.0, 1.0);
    let mut misses = Vec::new();
    for k in 1..200 {
        let f = k as f64 / 200.0;
        let split = surface.point_at(f * turns * tau, f * turns * pitch);
        match split_pcurve_at_image(&pc, &surface, &split, 1.0e-6) {
            Some((a, _)) => {
                let cut = a.point_at(a.t_end);
                if (cut.y - f * turns * pitch).abs() > 1.0e-6 {
                    misses.push(format!("f {f}: cut v {} on the wrong turn", cut.y));
                }
            }
            None => misses.push(format!("f {f}: split declined")),
        }
    }
    assert!(misses.is_empty(), "{} misses: {misses:?}", misses.len());
}
