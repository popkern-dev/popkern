//! PK-670 pins: the three graph-side closest-parameter searches on a
//! many-turn helical edge (radius / pitch ≈ 2.2, the M2.5 thread's
//! ratio) find the TRUE turn at 200 positions. Each failed on the
//! fixed-count code it replaced (negative controls recorded on the
//! PK-670 ledger), which aliased onto adjacent turns.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::fitting::fit_curve;
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{CurveOnEdge, Point3, Vec3};
use std::f64::consts::TAU;

const TURNS: f64 = 60.0;
const PITCH: f64 = 0.45;

fn cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    }
}

/// The helix as the cylinder's own chart image: u = θ, v = P·θ/2π.
fn helix_uv(theta: f64) -> (f64, f64) {
    (theta, PITCH * theta / TAU)
}

/// A NURBS helix of `TURNS` turns on [0, 1] (8 control points per
/// turn, fitted through 32 exact samples per turn — the swept thread
/// edges carry a comparable span density).
fn nurbs_helix() -> CurveOnEdge {
    let surface = cylinder();
    let m = (32.0 * TURNS) as usize + 1;
    let params: Vec<f64> = (0..m).map(|i| i as f64 / (m - 1) as f64).collect();
    let points: Vec<Point3> = params
        .iter()
        .map(|&s| {
            let (u, v) = helix_uv(s * TURNS * TAU);
            surface.point_at(u, v)
        })
        .collect();
    let n = (8.0 * TURNS) as usize;
    let curve = fit_curve(&points, &params, 3, n).expect("fit");
    CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve),
        t_start: 0.0,
        t_end: 1.0,
    }
}

#[test]
fn window_split_lands_on_the_true_turn() {
    let edge = nurbs_helix();
    let mut misses = Vec::new();
    for k in 0..200 {
        let t_true = (k as f64 + 0.5) / 200.0;
        let pos = edge.geometry.point_at(t_true);
        let t = super::super::graph_split_node::curve_window_interior_split(&edge, &pos)
            .expect("interior");
        if (t - t_true).abs() > 1e-6 {
            misses.push((t_true, t));
        }
    }
    assert!(
        misses.is_empty(),
        "{} of 200 splits off-turn: {:?}",
        misses.len(),
        &misses[..misses.len().min(5)]
    );
}

#[test]
fn projection_run_embeds_every_turn() {
    let edge = nurbs_helix();
    let surface = cylinder();
    let mut misses = Vec::new();
    for k in 0..200 {
        // 42-turn windows sliding along the helix.
        let a = 0.25 * k as f64 / 200.0;
        let window = CurveOnEdge {
            geometry: edge.geometry.clone(),
            t_start: a,
            t_end: a + 0.7,
        };
        let (run, _, _) =
            super::super::walk::chart_embedding::projection_run(&window, true, &surface)
                .expect("run");
        let du = run[run.len() - 1].0 - run[0].0;
        let want = 0.7 * TURNS * TAU;
        if (du - want).abs() > 1e-3 {
            misses.push((a, du / TAU));
        }
    }
    assert!(
        misses.is_empty(),
        "{} of 200 runs mis-wound (turns): {:?}",
        misses.len(),
        &misses[..misses.len().min(5)]
    );
}

#[test]
fn pcurve_station_finds_the_true_turn() {
    let surface = cylinder();
    let (u1, v1) = helix_uv(TURNS * TAU);
    let pc = Pcurve {
        geometry: PcurveGeometry::Linear {
            start: (0.0, 0.0),
            end: (u1, v1),
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    let mut misses = Vec::new();
    for k in 0..200 {
        let t_true = (k as f64 + 0.37) / 200.0;
        let (u, v) = helix_uv(t_true * TURNS * TAU);
        let p = surface.point_at(u, v);
        let (d, t) = super::super::graph_split_pcurve::pcurve_station_near(&pc, &surface, &p, 1.0);
        if d > 1e-3 || (t - t_true).abs() > 1e-4 {
            misses.push((t_true, t, d));
        }
    }
    assert!(
        misses.is_empty(),
        "{} of 200 stations off-turn: {:?}",
        misses.len(),
        &misses[..misses.len().min(5)]
    );
}
