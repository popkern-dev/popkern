//! Diagnostic support-frame rewrite for injected Gregory cap boundaries.
//!
//! The caller keeps this env-gated. The rewrite derives seam ownership
//! from the actual cap face loop and cap-side pcurves, then changes
//! only first-interior control points so boundary curves remain fixed.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

pub(super) fn apply_face_owned_support_frames(brep: &mut BRep, cap_face: FaceHandle) {
    let support_sides = collect_support_sides(brep, cap_face);
    if support_sides.is_empty() {
        return;
    }
    let trace = crate::a2_trace::support_frames_enabled();
    let Some(surface) = brep
        .faces
        .get_mut(cap_face)
        .and_then(|f| f.surface.as_mut())
    else {
        return;
    };
    let SurfaceGeometry::Gregory(patch) = &mut surface.geometry else {
        return;
    };
    let support_side_set = SupportSideSet::from_support_sides(&support_sides);
    for support_side in support_sides {
        if trace {
            eprintln!(
                "[a2 support-frame] face={cap_face:?} variant={} side={}",
                patch_variant_name(patch),
                support_side.side.name()
            );
        }
        match patch {
            GregoryPatch::RectangularBiquinticFull(p) => {
                apply_side_degree5(
                    &mut p.grid,
                    support_side.side,
                    &support_side.support,
                    &support_side.twin_pcurve,
                    support_side_set,
                );
            }
            GregoryPatch::RectangularBiseptic(p) => {
                apply_side_degree7(
                    &mut p.grid,
                    support_side.side,
                    &support_side.support,
                    &support_side.twin_pcurve,
                    support_side_set,
                );
            }
            GregoryPatch::RectangularTwoStripe(p) => {
                apply_side_degree7(
                    &mut p.grid,
                    support_side.side,
                    &support_side.support,
                    &support_side.twin_pcurve,
                    support_side_set,
                );
            }
            GregoryPatch::RectangularBiNonic(p) => {
                apply_side_degree9(
                    &mut p.grid,
                    support_side.side,
                    &support_side.support,
                    &support_side.twin_pcurve,
                    support_side_set,
                );
            }
            _ => {}
        }
    }
}

fn patch_variant_name(patch: &GregoryPatch) -> &'static str {
    match patch {
        GregoryPatch::Rectangular(_) => "rectangular",
        GregoryPatch::RectangularDeg5(_) => "rectangular_deg5",
        GregoryPatch::RectangularBiquinticFull(_) => "biquintic_full",
        GregoryPatch::RectangularBiseptic(_) => "biseptic",
        GregoryPatch::RectangularTwoStripe(_) => "two_stripe",
        GregoryPatch::RectangularBiNonic(_) => "binonic",
        GregoryPatch::Triangular(_) => "triangular",
    }
}

struct SupportSide {
    side: BoundarySide,
    support: SurfaceGeometry,
    twin_pcurve: Pcurve,
}

#[derive(Clone, Copy)]
enum BoundarySide {
    UZero,
    UOne,
    VZero,
    VOne,
}

impl BoundarySide {
    fn name(self) -> &'static str {
        match self {
            BoundarySide::UZero => "u0",
            BoundarySide::UOne => "u1",
            BoundarySide::VZero => "v0",
            BoundarySide::VOne => "v1",
        }
    }
}

#[derive(Clone, Copy, Default)]
struct SupportSideSet {
    u_zero: bool,
    u_one: bool,
    v_zero: bool,
    v_one: bool,
}

impl SupportSideSet {
    fn from_support_sides(sides: &[SupportSide]) -> Self {
        let mut out = Self::default();
        for side in sides {
            match side.side {
                BoundarySide::UZero => out.u_zero = true,
                BoundarySide::UOne => out.u_one = true,
                BoundarySide::VZero => out.v_zero = true,
                BoundarySide::VOne => out.v_one = true,
            }
        }
        out
    }

    fn contains(self, side: BoundarySide) -> bool {
        match side {
            BoundarySide::UZero => self.u_zero,
            BoundarySide::UOne => self.u_one,
            BoundarySide::VZero => self.v_zero,
            BoundarySide::VOne => self.v_one,
        }
    }
}

fn collect_support_sides(brep: &BRep, cap_face: FaceHandle) -> Vec<SupportSide> {
    let Some(face_data) = brep.faces.get(cap_face) else {
        return Vec::new();
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return Vec::new();
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut sides = Vec::new();
    for _ in 0..(brep.half_edges.len() + 1) {
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        if let Some(pcurve) = he.pcurve.as_ref()
            && let Some(side) = classify_boundary_side(pcurve)
        {
            let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
            if let Some(support) = brep
                .faces
                .get(twin_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.clone())
                && let Some(twin_pcurve) = brep.half_edges[he.twin].pcurve.clone()
            {
                sides.push(SupportSide {
                    side,
                    support,
                    twin_pcurve,
                });
            }
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    sides
}

fn classify_boundary_side(pcurve: &Pcurve) -> Option<BoundarySide> {
    let (t0, t1) = pcurve.parameter_range();
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;
    for i in 0..5 {
        let s = i as f64 / 4.0;
        let t = t0 + s * (t1 - t0);
        let uv = pcurve.point_at(t);
        u_min = u_min.min(uv.x);
        u_max = u_max.max(uv.x);
        v_min = v_min.min(uv.y);
        v_max = v_max.max(uv.y);
    }
    let tol = 1e-6;
    if u_min.abs() <= tol && u_max.abs() <= tol {
        Some(BoundarySide::UZero)
    } else if (u_min - 1.0).abs() <= tol && (u_max - 1.0).abs() <= tol {
        Some(BoundarySide::UOne)
    } else if v_min.abs() <= tol && v_max.abs() <= tol {
        Some(BoundarySide::VZero)
    } else if (v_min - 1.0).abs() <= tol && (v_max - 1.0).abs() <= tol {
        Some(BoundarySide::VOne)
    } else {
        None
    }
}

fn apply_side_degree7(
    grid: &mut [[Point3; 8]; 8],
    side: BoundarySide,
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    support_side_set: SupportSideSet,
) {
    apply_side(
        grid,
        side,
        support,
        twin_pcurve,
        SUPPORT_FRAME_STEP_SCALE / 7.0,
        support_side_set,
    );
}

fn apply_side_degree5(
    grid: &mut [[Point3; 6]; 6],
    side: BoundarySide,
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    support_side_set: SupportSideSet,
) {
    apply_side(
        grid,
        side,
        support,
        twin_pcurve,
        SUPPORT_FRAME_STEP_SCALE / 5.0,
        support_side_set,
    );
}

fn apply_side_degree9(
    grid: &mut [[Point3; 10]; 10],
    side: BoundarySide,
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    support_side_set: SupportSideSet,
) {
    apply_side(
        grid,
        side,
        support,
        twin_pcurve,
        SUPPORT_FRAME_STEP_SCALE / 9.0,
        support_side_set,
    );
}

fn apply_side<const N: usize>(
    grid: &mut [[Point3; N]; N],
    side: BoundarySide,
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    step: f64,
    support_side_set: SupportSideSet,
) {
    let mut boundary = [Point3::new(0.0, 0.0, 0.0); N];
    match side {
        BoundarySide::UZero => {
            boundary.copy_from_slice(&grid[0]);
            let cross_cps = selected_cross_tangent_cps(&boundary, support, twin_pcurve, side);
            let j_start = if support_side_set.contains(BoundarySide::VZero) {
                2
            } else {
                0
            };
            let j_end = if support_side_set.contains(BoundarySide::VOne) {
                N - 2
            } else {
                N
            };
            for j in j_start..j_end {
                grid[1][j] = offset_point(boundary[j], cross_cps[j], step);
            }
        }
        BoundarySide::UOne => {
            boundary.copy_from_slice(&grid[N - 1]);
            let cross_cps = selected_cross_tangent_cps(&boundary, support, twin_pcurve, side);
            let j_start = if support_side_set.contains(BoundarySide::VZero) {
                2
            } else {
                0
            };
            let j_end = if support_side_set.contains(BoundarySide::VOne) {
                N - 2
            } else {
                N
            };
            for j in j_start..j_end {
                grid[N - 2][j] = offset_point(boundary[j], cross_cps[j], -step);
            }
        }
        BoundarySide::VZero => {
            for i in 0..N {
                boundary[i] = grid[i][0];
            }
            let cross_cps = selected_cross_tangent_cps(&boundary, support, twin_pcurve, side);
            let i_start = if support_side_set.contains(BoundarySide::UZero) {
                2
            } else {
                0
            };
            let i_end = if support_side_set.contains(BoundarySide::UOne) {
                N - 2
            } else {
                N
            };
            for i in i_start..i_end {
                grid[i][1] = offset_point(boundary[i], cross_cps[i], step);
            }
        }
        BoundarySide::VOne => {
            for i in 0..N {
                boundary[i] = grid[i][N - 1];
            }
            let cross_cps = selected_cross_tangent_cps(&boundary, support, twin_pcurve, side);
            let i_start = if support_side_set.contains(BoundarySide::UZero) {
                2
            } else {
                0
            };
            let i_end = if support_side_set.contains(BoundarySide::UOne) {
                N - 2
            } else {
                N
            };
            for i in i_start..i_end {
                grid[i][N - 2] = offset_point(boundary[i], cross_cps[i], -step);
            }
        }
    }
}

const SUPPORT_FRAME_STEP_SCALE: f64 = 1.0;

fn selected_cross_tangent_cps<const N: usize>(
    boundary: &[Point3; N],
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    side: BoundarySide,
) -> [Vec3; N] {
    let fitted = fitted_cross_tangent_cps(boundary, support, twin_pcurve, side);
    let local = cp_local_cross_tangent_cps(boundary, support, twin_pcurve, side);
    let fitted_score = score_cross_tangent_cps(boundary, support, twin_pcurve, side, &fitted);
    let local_score = score_cross_tangent_cps(boundary, support, twin_pcurve, side, &local);
    if crate::a2_trace::support_frames_enabled() {
        let trace_t = 0.65;
        let trace_point = bezier_point(boundary, trace_t);
        let trace_tangent = bezier_tangent(boundary, trace_t);
        let trace_normal =
            support_normal_at_boundary_t(support, twin_pcurve, &trace_point, trace_t);
        let trace_cross = bezier_vec(&fitted, trace_t);
        eprintln!(
            "[a2 support-frame] side={} fitted_score={fitted_score:.6e} local_score={local_score:.6e} t=.65 tangent=({:.3},{:.3},{:.3}) normal={} fitted_cross=({:.3},{:.3},{:.3})",
            side.name(),
            trace_tangent.x,
            trace_tangent.y,
            trace_tangent.z,
            trace_normal
                .map(|n| format!("({:.3},{:.3},{:.3})", n.x, n.y, n.z))
                .unwrap_or_else(|| "none".to_string()),
            trace_cross.x,
            trace_cross.y,
            trace_cross.z
        );
    }
    if fitted_score <= local_score {
        fitted
    } else {
        local
    }
}

fn fitted_cross_tangent_cps<const N: usize>(
    boundary: &[Point3; N],
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    side: BoundarySide,
) -> [Vec3; N] {
    use popkern_core::geometry::fitting::fit_curve;

    let degree = N - 1;
    let n_samples = 40usize.max(N * 4);
    let mut samples = Vec::with_capacity(n_samples);
    let mut params = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1) as f64;
        let point = bezier_point(boundary, t);
        let tangent = bezier_tangent(boundary, t);
        let Some(normal) = support_normal_at_boundary_t(support, twin_pcurve, &point, t) else {
            return cp_local_cross_tangent_cps(boundary, support, twin_pcurve, side);
        };
        let cross = oriented_cross_tangent(tangent, normal, side)
            .normalized()
            .map(|unit| unit * tangent.length())
            .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
        samples.push(Point3::new(cross.x, cross.y, cross.z));
        params.push(t);
    }
    let Ok(fitted) = fit_curve(&samples, &params, degree, N) else {
        return cp_local_cross_tangent_cps(boundary, support, twin_pcurve, side);
    };
    if fitted.control_points.len() != N {
        return cp_local_cross_tangent_cps(boundary, support, twin_pcurve, side);
    }
    let mut out = [Vec3::new(0.0, 0.0, 0.0); N];
    for (i, cp) in fitted.control_points.iter().enumerate() {
        out[i] = Vec3::new(cp.x, cp.y, cp.z);
    }
    out
}

fn score_cross_tangent_cps<const N: usize>(
    boundary: &[Point3; N],
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    side: BoundarySide,
    cross_cps: &[Vec3; N],
) -> f64 {
    let n_samples = 24usize;
    let mut worst = 0.0_f64;
    for i in 0..n_samples {
        let t = 0.025 + 0.95 * (i as f64) / (n_samples - 1) as f64;
        let point = bezier_point(boundary, t);
        let tangent = bezier_tangent(boundary, t);
        let cross = bezier_vec(cross_cps, t);
        let Some(support_normal) = support_normal_at_boundary_t(support, twin_pcurve, &point, t)
        else {
            return f64::INFINITY;
        };
        let cap_normal = match side {
            BoundarySide::UZero | BoundarySide::UOne => cross.cross(&tangent),
            BoundarySide::VZero | BoundarySide::VOne => tangent.cross(&cross),
        };
        let Some(cap_unit) = cap_normal.normalized() else {
            return f64::INFINITY;
        };
        let Some(support_unit) = support_normal.normalized() else {
            return f64::INFINITY;
        };
        worst = worst.max(1.0 - cap_unit.dot(&support_unit).abs());
    }
    worst
}

fn cp_local_cross_tangent_cps<const N: usize>(
    boundary: &[Point3; N],
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    side: BoundarySide,
) -> [Vec3; N] {
    let mut out = [Vec3::new(0.0, 0.0, 0.0); N];
    for (i, item) in out.iter_mut().enumerate() {
        let t = if N > 1 {
            i as f64 / (N - 1) as f64
        } else {
            0.0
        };
        if let Some(cross) =
            boundary_owned_cross_tangent(boundary, i, support, twin_pcurve, t, side)
        {
            *item = cross;
        }
    }
    out
}

fn bezier_vec<const N: usize>(cps: &[Vec3; N], t: f64) -> Vec3 {
    let degree = N - 1;
    let mut out = Vec3::new(0.0, 0.0, 0.0);
    for (i, cp) in cps.iter().enumerate() {
        out = out + *cp * bernstein(degree, i, t);
    }
    out
}

fn bezier_point<const N: usize>(cps: &[Point3; N], t: f64) -> Point3 {
    let degree = N - 1;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for (i, cp) in cps.iter().enumerate() {
        let b = bernstein(degree, i, t);
        x += b * cp.x;
        y += b * cp.y;
        z += b * cp.z;
    }
    Point3::new(x, y, z)
}

fn bezier_tangent<const N: usize>(cps: &[Point3; N], t: f64) -> Vec3 {
    if N < 2 {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    let degree = N - 1;
    let deriv_degree = degree - 1;
    let mut out = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..degree {
        let delta = point_delta(cps[i + 1], cps[i]) * degree as f64;
        out = out + delta * bernstein(deriv_degree, i, t);
    }
    out
}

fn bernstein(degree: usize, i: usize, t: f64) -> f64 {
    binomial(degree, i) * t.powi(i as i32) * (1.0 - t).powi((degree - i) as i32)
}

fn binomial(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    let mut out = 1.0;
    for i in 0..k {
        out *= (n - i) as f64 / (i + 1) as f64;
    }
    out
}

fn boundary_owned_cross_tangent<const N: usize>(
    boundary: &[Point3; N],
    index: usize,
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    t: f64,
    side: BoundarySide,
) -> Option<Vec3> {
    let tangent = boundary_cp_tangent(boundary, index)?;
    let normal = support_normal_at_boundary_t(support, twin_pcurve, &boundary[index], t)?;
    let cross = oriented_cross_tangent(tangent, normal, side);
    cross.normalized().map(|unit| unit * tangent.length())
}

fn oriented_cross_tangent(tangent: Vec3, normal: Vec3, side: BoundarySide) -> Vec3 {
    match side {
        BoundarySide::UZero | BoundarySide::UOne => tangent.cross(&normal),
        BoundarySide::VZero | BoundarySide::VOne => normal.cross(&tangent),
    }
}

fn boundary_cp_tangent<const N: usize>(boundary: &[Point3; N], index: usize) -> Option<Vec3> {
    if N < 2 || index >= N {
        return None;
    }
    let degree = (N - 1) as f64;
    let tangent = if index == 0 {
        point_delta(boundary[1], boundary[0]) * degree
    } else if index + 1 == N {
        point_delta(boundary[N - 1], boundary[N - 2]) * degree
    } else {
        point_delta(boundary[index + 1], boundary[index - 1]) * (0.5 * degree)
    };
    (tangent.length() > 1e-12).then_some(tangent)
}

fn support_normal_at_boundary_t(
    support: &SurfaceGeometry,
    twin_pcurve: &Pcurve,
    point: &Point3,
    boundary_t: f64,
) -> Option<Vec3> {
    let (t0, t1) = twin_pcurve.parameter_range();
    let forward_t = t0 + boundary_t * (t1 - t0);
    let reverse_t = t1 - boundary_t * (t1 - t0);
    let forward_uv = twin_pcurve.point_at(forward_t);
    let reverse_uv = twin_pcurve.point_at(reverse_t);
    let forward_p = support.point_at(forward_uv.x, forward_uv.y);
    let reverse_p = support.point_at(reverse_uv.x, reverse_uv.y);
    let uv = if point.distance_to(&forward_p) <= point.distance_to(&reverse_p) {
        forward_uv
    } else {
        reverse_uv
    };
    support.normal_at(uv.x, uv.y).normalized()
}

fn point_delta(a: Point3, b: Point3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

fn offset_point(point: Point3, direction: Vec3, scale: f64) -> Point3 {
    Point3::new(
        point.x + scale * direction.x,
        point.y + scale * direction.y,
        point.z + scale * direction.z,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plane_z() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    fn linear_pcurve() -> Pcurve {
        Pcurve::linear((0.0, 0.0), (1.0, 1.0), 0.0, 1.0)
    }

    fn regular_grid_10() -> [[Point3; 10]; 10] {
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
        for (i, row) in grid.iter_mut().enumerate() {
            for (j, cp) in row.iter_mut().enumerate() {
                *cp = Point3::new(i as f64, j as f64, 0.0);
            }
        }
        grid
    }

    fn assert_point_close(actual: Point3, expected: Point3) {
        assert!(
            actual.distance_to(&expected) <= 1e-9,
            "actual={actual:?} expected={expected:?}"
        );
    }

    #[test]
    fn support_frame_preserves_boundary_rows_and_columns() {
        let mut grid = regular_grid_10();
        let before = grid;
        let side_set = SupportSideSet {
            u_zero: true,
            u_one: false,
            v_zero: false,
            v_one: false,
        };

        apply_side(
            &mut grid,
            BoundarySide::UZero,
            &plane_z(),
            &linear_pcurve(),
            1.0 / 9.0,
            side_set,
        );

        assert_eq!(grid[0], before[0], "u=0 boundary row moved");
        for i in 0..10 {
            assert_eq!(grid[i][0], before[i][0], "v=0 boundary column moved");
            assert_eq!(grid[i][9], before[i][9], "v=1 boundary column moved");
        }
        assert_point_close(grid[1][5], Point3::new(1.0, 5.0, 0.0));
    }

    #[test]
    fn adjacent_support_sides_do_not_overwrite_each_other_first_strips() {
        let mut grid = regular_grid_10();
        let side_set = SupportSideSet {
            u_zero: true,
            u_one: false,
            v_zero: false,
            v_one: true,
        };

        apply_side(
            &mut grid,
            BoundarySide::UZero,
            &plane_z(),
            &linear_pcurve(),
            1.0 / 9.0,
            side_set,
        );
        let u_zero_owned = grid[1][7];
        let shared_corner_strip = grid[1][8];

        apply_side(
            &mut grid,
            BoundarySide::VOne,
            &plane_z(),
            &linear_pcurve(),
            1.0 / 9.0,
            side_set,
        );

        assert_point_close(grid[1][7], u_zero_owned);
        assert_point_close(grid[1][8], shared_corner_strip);
        assert_point_close(grid[5][8], Point3::new(5.0, 8.0, 0.0));
    }

    #[test]
    fn support_normals_are_read_from_twin_pcurve_parameterization() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let twin_pcurve = Pcurve::linear((0.0, 0.0), (0.0, 1.0), 0.0, 1.0);
        let misleading_point = cyl.point_at(std::f64::consts::PI, 0.5);

        let normal = support_normal_at_boundary_t(&cyl, &twin_pcurve, &misleading_point, 0.5)
            .expect("normal from twin pcurve");
        let expected = cyl.normal_at(0.0, 0.5).normalized().unwrap();

        assert!(
            (normal - expected).length() <= 1e-9,
            "normal={normal:?} expected={expected:?}"
        );
    }
}
