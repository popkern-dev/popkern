//! Diagnostic support-frame rewrite for injected Gregory cap boundaries.
//!
//! The caller keeps this env-gated. The rewrite derives seam ownership
//! from the actual cap face loop and cap-side pcurves, then changes
//! only first-interior control points so boundary curves remain fixed.

#[path = "boundary_support_cross_tangent.rs"]
mod boundary_support_cross_tangent;

use boundary_support_cross_tangent::{offset_point, selected_cross_tangent_cps};
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

pub(super) fn apply_face_owned_support_frames(brep: &mut BRep, cap_face: FaceHandle) {
    let support_sides = collect_support_sides(brep, cap_face);
    if support_sides.is_empty() {
        return;
    }
    let trace = crate::direct_apex_trace::support_frames_enabled();
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
                "[direct-apex support-frame] face={cap_face:?} variant={} side={}",
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
        GregoryPatch::RectangularBiNonicDualTwist(_) => "binonic_dual_twist",
        GregoryPatch::Triangular(_) => "triangular",
    }
}

struct SupportSide {
    side: BoundarySide,
    support: SurfaceGeometry,
    twin_pcurve: Pcurve,
}

#[derive(Clone, Copy)]
pub(super) enum BoundarySide {
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

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::surface::SurfaceEval;

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

        let normal = boundary_support_cross_tangent::support_normal_at_boundary_t(
            &cyl,
            &twin_pcurve,
            &misleading_point,
            0.5,
        )
        .expect("normal from twin pcurve");
        let expected = cyl.normal_at(0.0, 0.5).normalized().unwrap();

        assert!(
            (normal - expected).length() <= 1e-9,
            "normal={normal:?} expected={expected:?}"
        );
    }
}
