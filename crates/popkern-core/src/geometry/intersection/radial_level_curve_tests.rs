//! Unit tests for the radial level-curve fast path (PK-461).
//!
//! Fixture: an exact ruled helical band between two stitched
//! quarter-arc rational helices at radii r1 < r2 sharing knots and
//! weights. On that surface the radial field is exactly
//! `ρ(u, v) = (1−u)·r1 + u·r2` (the outer row is a radial scaling of
//! the inner and the shared weights make the ruling a Euclidean
//! lerp), so the true intersection with a coaxial cylinder at
//! r1 < R < r2 is the v-isocurve at `u* = (R−r1)/(r2−r1)` — a known
//! curve the emitted records must track.

use std::f64::consts::{FRAC_1_SQRT_2, TAU};
use std::sync::{Mutex, MutexGuard, OnceLock};

use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::radial_level_curve::try_radial_level_curve;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// The kill-switch test mutates process environment; every test in
/// this file serializes on this lock so the arm-must-fire assertions
/// never race the pinned-marching window.
pub(super) fn env_lock() -> MutexGuard<'static, ()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    match LOCK.get_or_init(|| Mutex::new(())).lock() {
        Ok(g) => g,
        Err(poisoned) => poisoned.into_inner(),
    }
}

/// Quarter-arc rational helix control row at `radius` (the
/// popkern-sweep canonical constant helix spelling).
fn helix_row(radius: f64, turns: f64, pitch: f64, start_angle: f64) -> (Vec<Point3>, Vec<f64>) {
    let n_quarters = (turns * 4.0).ceil() as usize;
    let total_angle = turns * TAU;
    let total_height = turns * pitch;
    let mut points = Vec::new();
    let mut weights = Vec::new();
    for i in 0..=n_quarters {
        let frac = i as f64 / n_quarters as f64;
        let angle = start_angle + frac * total_angle;
        let z = frac * total_height;
        points.push(Point3::new(radius * angle.cos(), radius * angle.sin(), z));
        weights.push(1.0);
        if i < n_quarters {
            let frac_mid = (i as f64 + 0.5) / n_quarters as f64;
            let angle_mid = start_angle + frac_mid * total_angle;
            let z_mid = frac_mid * total_height;
            let scale = radius / FRAC_1_SQRT_2;
            points.push(Point3::new(
                scale * angle_mid.cos(),
                scale * angle_mid.sin(),
                z_mid,
            ));
            weights.push(FRAC_1_SQRT_2);
        }
    }
    (points, weights)
}

fn quarter_knots(turns: f64) -> Vec<f64> {
    let n_quarters = (turns * 4.0).ceil() as usize;
    let mut knots = vec![0.0, 0.0, 0.0];
    for i in 1..n_quarters {
        let k = i as f64 / n_quarters as f64;
        knots.push(k);
        knots.push(k);
    }
    knots.extend([1.0, 1.0, 1.0]);
    knots
}

/// Ruled band: u ∈ [0,1] radial (r1 → r2), v along the helix.
pub(super) fn helical_band(r1: f64, r2: f64, turns: f64, pitch: f64) -> SurfaceGeometry {
    let (inner, w_inner) = helix_row(r1, turns, pitch, 0.0);
    let (outer, _) = helix_row(r2, turns, pitch, 0.0);
    let num_v = inner.len();
    // Grid indexing is control_points[i_u * num_v + i_v].
    let mut control_points = inner;
    control_points.extend(outer);
    let mut weights = w_inner.clone();
    weights.extend(w_inner);
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: quarter_knots(turns),
        control_points,
        num_u: 2,
        num_v,
        weights,
    })
}

fn shank_cylinder(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

pub(super) fn curves(
    result: Vec<SurfaceSurfaceIntersection>,
) -> Vec<crate::geometry::curve_on_surface::CurveOnSurface> {
    result
        .into_iter()
        .map(|i| match i {
            SurfaceSurfaceIntersection::Curve(c) => c,
            other => panic!("expected curve record, got {other:?}"),
        })
        .collect()
}

#[test]
fn helical_band_crossing_emits_per_wrap_records() {
    let _guard = env_lock();
    let turns = 3.0;
    let pitch = 1.0;
    let (r1, r2, level) = (2.0, 4.0, 3.0);
    let wall = helical_band(r1, r2, turns, pitch);
    let cyl = shank_cylinder(level);
    let tol = Tolerance::default();
    let result = try_radial_level_curve(
        &wall,
        (0.0, 1.0),
        (0.0, 1.0),
        &cyl,
        (0.0, TAU),
        (-1.0, 4.0),
        &tol,
    )
    .expect("arm claims the multi-wrap band pair");
    let records = curves(result);
    // Three wraps split at the chart seam into three records (the
    // trace starts and ends exactly on the seam; end blips are
    // dropped by the chord filter).
    assert_eq!(records.len(), 3, "one record per wrap");

    let band = 10.0 * tol.linear();
    let mut total_len = 0.0;
    for record in &records {
        let (t0, t1) = record.curve_3d.parameter_range();
        let mut prev: Option<Point3> = None;
        for i in 0..=64 {
            let t = t0 + (t1 - t0) * i as f64 / 64.0;
            let p = record.curve_3d.point_at(t);
            // On the cylinder: |ρ − R| within the honesty band.
            let rho = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (rho - level).abs() <= band,
                "record point off the cylinder by {:.3e}",
                (rho - level).abs()
            );
            // Cylinder pcurve stays on the principal chart branch.
            let cuv = record.curve_uv2.point_at(t);
            assert!(
                (-1e-9..=TAU + 1e-9).contains(&cuv.x),
                "cylinder pcurve angle {} left the principal branch",
                cuv.x
            );
            if let Some(q) = prev {
                total_len += q.distance_to(&p);
            }
            prev = Some(p);
        }
        // Interior seam endpoints land exactly on the seam.
        let u2_start = record.curve_uv2.point_at(t0).x;
        let u2_end = record.curve_uv2.point_at(t1).x;
        for u2 in [u2_start, u2_end] {
            assert!(
                u2.abs() < 1e-6 || (u2 - TAU).abs() < 1e-6,
                "record endpoint angle {u2} is not on the chart seam"
            );
        }
    }
    // The three records together track the full helix at the level
    // radius: length within 1% of the analytic helix length.
    let helix_len = turns * ((TAU * level).powi(2) + pitch * pitch).sqrt();
    assert!(
        (total_len - helix_len).abs() < 0.01 * helix_len,
        "total record length {total_len:.4} vs analytic helix {helix_len:.4}"
    );
}

#[test]
fn operand_order_swaps_pcurve_slots() {
    let _guard = env_lock();
    let wall = helical_band(2.0, 4.0, 3.0, 1.0);
    let cyl = shank_cylinder(3.0);
    let tol = Tolerance::default();
    let result = try_radial_level_curve(
        &cyl,
        (0.0, TAU),
        (-1.0, 4.0),
        &wall,
        (0.0, 1.0),
        (0.0, 1.0),
        &tol,
    )
    .expect("arm claims the pair in either operand order");
    let records = curves(result);
    assert_eq!(records.len(), 3);
    for record in &records {
        let (t0, t1) = record.curve_3d.parameter_range();
        let mid = 0.5 * (t0 + t1);
        // Surface 1 is the cylinder now: its pcurve angle is on the
        // principal branch, and the wall pcurve sits in [0,1]².
        let cyl_uv = record.curve_uv1.point_at(mid);
        assert!((-1e-9..=TAU + 1e-9).contains(&cyl_uv.x));
        let wall_uv = record.curve_uv2.point_at(mid);
        assert!((-1e-6..=1.0 + 1e-6).contains(&wall_uv.x));
        assert!((-1e-6..=1.0 + 1e-6).contains(&wall_uv.y));
    }
}

#[test]
fn sub_wrap_trace_declines_to_marching() {
    let _guard = env_lock();
    // One turn wraps exactly 1.0 ≤ the 1.5-period scope gate.
    let wall = helical_band(2.0, 4.0, 1.0, 1.0);
    let cyl = shank_cylinder(3.0);
    let tol = Tolerance::default();
    assert!(
        try_radial_level_curve(
            &wall,
            (0.0, 1.0),
            (0.0, 1.0),
            &cyl,
            (0.0, TAU),
            (-1.0, 2.0),
            &tol,
        )
        .is_none(),
        "sub-wrap pairs stay on the marching path"
    );
}

#[test]
fn tangent_level_radius_declines() {
    let _guard = env_lock();
    // Level radius at the band's outer boundary: no transversal
    // crossing exists; the graze/structure gates must decline.
    let wall = helical_band(2.0, 4.0, 3.0, 1.0);
    let cyl = shank_cylinder(4.0);
    let tol = Tolerance::default();
    assert!(
        try_radial_level_curve(
            &wall,
            (0.0, 1.0),
            (0.0, 1.0),
            &cyl,
            (0.0, TAU),
            (-1.0, 4.0),
            &tol,
        )
        .is_none(),
        "tangent contact is not a transversal graph"
    );
}

#[test]
fn radius_beyond_control_hull_answers_empty_without_marching() {
    let _guard = env_lock();
    let wall = helical_band(2.0, 4.0, 3.0, 1.0);
    // The mid-arc control points reach r·√2; beyond that the convex
    // hull bound PROVES the wall never reaches the cylinder (PK-513):
    // the arm answers "no intersection" instead of handing the pair
    // to a seed hunt over the multi-turn chart.
    let cyl = shank_cylinder(8.0);
    let tol = Tolerance::default();
    let answer = try_radial_level_curve(
        &wall,
        (0.0, 1.0),
        (0.0, 1.0),
        &cyl,
        (0.0, TAU),
        (-1.0, 4.0),
        &tol,
    );
    assert!(answer.is_some_and(|curves| curves.is_empty()));
}

#[test]
fn kill_switch_pins_the_marching_path() {
    let _guard = env_lock();
    let wall = helical_band(2.0, 4.0, 3.0, 1.0);
    let cyl = shank_cylinder(3.0);
    let tol = Tolerance::default();
    // SAFETY-adjacent note: env mutation is serialized by env_lock()
    // and scoped to this assertion window (torus-arm precedent).
    unsafe {
        std::env::set_var("POPKERN_ANALYTIC_HELICAL_SSI", "0");
    }
    let pinned = try_radial_level_curve(
        &wall,
        (0.0, 1.0),
        (0.0, 1.0),
        &cyl,
        (0.0, TAU),
        (-1.0, 4.0),
        &tol,
    );
    unsafe {
        std::env::remove_var("POPKERN_ANALYTIC_HELICAL_SSI");
    }
    assert!(pinned.is_none(), "kill switch must pin marching");
    assert!(
        try_radial_level_curve(
            &wall,
            (0.0, 1.0),
            (0.0, 1.0),
            &cyl,
            (0.0, TAU),
            (-1.0, 4.0),
            &tol,
        )
        .is_some(),
        "arm resumes once the switch is cleared"
    );
}

#[test]
fn partial_crossing_structure_declines() {
    let _guard = env_lock();
    // Outer row tapers from below to above the level radius: the
    // crossing exists only over part of the stations, so the count
    // consistency gate declines.
    let turns = 3.0;
    let pitch = 1.0;
    let (inner, w_inner) = helix_row(2.0, turns, pitch, 0.0);
    let (outer_base, _) = helix_row(1.0, turns, pitch, 0.0);
    let n = outer_base.len();
    let outer: Vec<Point3> = outer_base
        .iter()
        .enumerate()
        .map(|(j, p)| {
            // Radial scale grows 2.4 → 3.6 along v (level at 3.0).
            let s = 2.4 + 1.2 * j as f64 / (n - 1) as f64;
            Point3::new(p.x * s, p.y * s, p.z)
        })
        .collect();
    let mut control_points = inner;
    control_points.extend(outer);
    let mut weights = w_inner.clone();
    weights.extend(w_inner);
    let wall = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: quarter_knots(turns),
        control_points,
        num_u: 2,
        num_v: n,
        weights,
    });
    let cyl = shank_cylinder(3.0);
    let tol = Tolerance::default();
    assert!(
        try_radial_level_curve(
            &wall,
            (0.0, 1.0),
            (0.0, 1.0),
            &cyl,
            (0.0, TAU),
            (-1.0, 4.0),
            &tol,
        )
        .is_none(),
        "inconsistent crossing structure declines to marching"
    );
}

#[test]
fn left_handed_band_serves_symmetrically() {
    let _guard = env_lock();
    // A left-hand helix (negative total angle) — build by mirroring
    // the y coordinate, which reverses the winding sense.
    let wall = {
        let SurfaceGeometry::Nurbs(mut ns) = helical_band(2.0, 4.0, 3.0, 1.0) else {
            panic!("band is NURBS");
        };
        for p in &mut ns.control_points {
            p.y = -p.y;
        }
        SurfaceGeometry::Nurbs(ns)
    };
    let cyl = shank_cylinder(3.0);
    let tol = Tolerance::default();
    let result = try_radial_level_curve(
        &wall,
        (0.0, 1.0),
        (0.0, 1.0),
        &cyl,
        (0.0, TAU),
        (-1.0, 4.0),
        &tol,
    )
    .expect("left-hand winding serves");
    assert_eq!(curves(result).len(), 3);
}

#[test]
fn axial_window_clips_the_trace() {
    let _guard = env_lock();
    // Cylinder window covers only the middle third of the coil
    // (plus the arm's own pad of half the window span): the emitted
    // records must not reach far beyond the padded window.
    let wall = helical_band(2.0, 4.0, 6.0, 1.0); // z ∈ [0, 6]
    let cyl = shank_cylinder(3.0);
    let tol = Tolerance::default();
    let window = (2.5, 3.5);
    let result = try_radial_level_curve(
        &wall,
        (0.0, 1.0),
        (0.0, 1.0),
        &cyl,
        (0.0, TAU),
        window,
        &tol,
    )
    .expect("clipped trace still serves");
    let records = curves(result);
    assert!(!records.is_empty());
    let pad = 0.5 * (window.1 - window.0) + tol.linear();
    for record in &records {
        let (t0, t1) = record.curve_3d.parameter_range();
        for i in 0..=32 {
            let t = t0 + (t1 - t0) * i as f64 / 32.0;
            let z = record.curve_3d.point_at(t).z;
            assert!(
                z > window.0 - pad - 0.2 && z < window.1 + pad + 0.2,
                "record reaches z = {z} far outside the padded window"
            );
        }
    }
}
