//! Segment-aware G0 reconciliation for KP20 support-cell shared full boundaries.

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot;

#[derive(Debug, Clone, Copy)]
struct FullBoundaryTrim {
    face_a: FaceHandle,
    side_a: BoundarySide,
    start_a: f64,
    end_a: f64,
    face_b: FaceHandle,
    side_b: BoundarySide,
    start_b: f64,
    end_b: f64,
}

#[derive(Debug, Clone, Copy)]
enum BoundarySide {
    UMin,
    UMax,
    VMin,
    VMax,
}

pub(crate) fn reconcile_support_cell_full_boundaries(
    brep: &mut BRep,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
) -> Result<(), CoreError> {
    let trims = full_boundary_trims(brep, slots)?;
    for trim in trims {
        reconcile_full_boundary_trim(brep, trim)?;
    }
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell boundary reconciliation produced invalid BRep: {error}"
        ))
    })
}

pub(crate) fn reconcile_assigned_support_cell_full_boundaries(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<(), CoreError> {
    let mut slots = Vec::with_capacity(faces.len());
    for &face in faces {
        let Ok(half_edges) = face_outer_half_edges(brep, face) else {
            continue;
        };
        let boundary_edge_count = half_edges.len();
        let assignment = if boundary_edge_count == 4 {
            crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceAssignment::FullDomainRectangularPatch
        } else {
            crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace
        };
        slots.push(Kp20ExpandedSupportSurfaceSlot {
            face,
            boundary_edge_count,
            assignment,
        });
    }
    reconcile_support_cell_full_boundaries(brep, &slots)
}

fn full_boundary_trims(
    brep: &BRep,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
) -> Result<Vec<FullBoundaryTrim>, CoreError> {
    let mut out = Vec::new();
    for (slot_a, slot) in slots.iter().enumerate() {
        for &half_edge in &face_outer_half_edges(brep, slot.face)? {
            let twin = brep.half_edges[half_edge].twin;
            let Some(face_b) = brep
                .loops
                .get(brep.half_edges[twin].loop_)
                .map(|loop_| loop_.face)
            else {
                continue;
            };
            let Some(slot_b) = slots.iter().position(|slot| slot.face == face_b) else {
                continue;
            };
            if slot_b <= slot_a {
                continue;
            }
            let Some(pcurve_a) = brep.half_edges[half_edge].pcurve.as_ref() else {
                continue;
            };
            let Some(pcurve_b) = brep.half_edges[twin].pcurve.as_ref() else {
                continue;
            };
            let Some((side_a, start_a, end_a)) = boundary_span(pcurve_a) else {
                continue;
            };
            let Some((side_b, start_b, end_b)) = boundary_span(pcurve_b) else {
                continue;
            };
            out.push(FullBoundaryTrim {
                face_a: slot.face,
                side_a,
                start_a,
                end_a,
                face_b,
                side_b,
                start_b,
                end_b,
            });
        }
    }
    Ok(out)
}

fn reconcile_full_boundary_trim(brep: &mut BRep, trim: FullBoundaryTrim) -> Result<(), CoreError> {
    let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
    let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
    let before_gap = trim_gap(&patch_a, &patch_b, trim);
    let mut candidate_a = patch_a.clone();
    let mut candidate_b = patch_b.clone();
    let mut params_a = [0.0; 10];
    let mut params_b = [0.0; 10];
    let mut targets = [Point3::new(0.0, 0.0, 0.0); 10];
    for index in 0..10 {
        let fraction = index as f64 / 9.0;
        params_a[index] = lerp_scalar(trim.start_a, trim.end_a, fraction);
        params_b[index] = lerp_scalar(trim.start_b, trim.end_b, 1.0 - fraction);
        targets[index] = average_point(
            eval_boundary(&patch_a, trim.side_a, params_a[index]),
            eval_boundary(&patch_b, trim.side_b, params_b[index]),
        );
    }
    fit_boundary_controls(&mut candidate_a, trim.side_a, &params_a, &targets);
    fit_boundary_controls(&mut candidate_b, trim.side_b, &params_b, &targets);
    let after_gap = trim_gap(&candidate_a, &candidate_b, trim);
    if after_gap + 1.0e-12 < before_gap
        && max_boundary_displacement(&patch_a, &candidate_a, trim.side_a) <= 0.1
        && max_boundary_displacement(&patch_b, &candidate_b, trim.side_b) <= 0.1
    {
        set_rectangular_binonic_on_face(brep, trim.face_a, candidate_a)?;
        set_rectangular_binonic_on_face(brep, trim.face_b, candidate_b)?;
    }
    Ok(())
}

fn rectangular_binonic_on_face(
    brep: &BRep,
    face: FaceHandle,
) -> Result<popkern_core::geometry::gregory::RectangularGregoryBiNonic, CoreError> {
    let Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(patch))) = brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Err(CoreError::PreconditionFailed(
            "kp20 support-cell boundary reconciliation requires rectangular bi-nonic surfaces"
                .into(),
        ));
    };
    Ok((**patch).clone())
}

fn set_rectangular_binonic_on_face(
    brep: &mut BRep,
    face: FaceHandle,
    patch: popkern_core::geometry::gregory::RectangularGregoryBiNonic,
) -> Result<(), CoreError> {
    let Some(face_data) = brep.faces.get_mut(face) else {
        return Err(CoreError::InvalidHandle(
            "kp20 support-cell boundary reconciliation face invalid".into(),
        ));
    };
    face_data.surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(Box::new(patch))),
        orientation: true,
    });
    Ok(())
}

fn boundary_span(pcurve: &Pcurve) -> Option<(BoundarySide, f64, f64)> {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    if close(start.x, 0.0) && close(end.x, 0.0) && usable_span(start.y, end.y) {
        Some((BoundarySide::UMin, start.y, end.y))
    } else if close(start.x, 1.0) && close(end.x, 1.0) && usable_span(start.y, end.y) {
        Some((BoundarySide::UMax, start.y, end.y))
    } else if close(start.y, 0.0) && close(end.y, 0.0) && usable_span(start.x, end.x) {
        Some((BoundarySide::VMin, start.x, end.x))
    } else if close(start.y, 1.0) && close(end.y, 1.0) && usable_span(start.x, end.x) {
        Some((BoundarySide::VMax, start.x, end.x))
    } else {
        None
    }
}

fn usable_span(a: f64, b: f64) -> bool {
    (a - b).abs() > 1.0e-9 && (0.0..=1.0).contains(&a) && (0.0..=1.0).contains(&b)
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}

fn eval_boundary(
    patch: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    side: BoundarySide,
    t: f64,
) -> Point3 {
    match side {
        BoundarySide::UMin => patch.point_at(0.0, t),
        BoundarySide::UMax => patch.point_at(1.0, t),
        BoundarySide::VMin => patch.point_at(t, 0.0),
        BoundarySide::VMax => patch.point_at(t, 1.0),
    }
}

fn fit_boundary_controls(
    patch: &mut popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    side: BoundarySide,
    params: &[f64; 10],
    targets: &[Point3; 10],
) {
    let prior = boundary_controls(patch, side);
    let lambda = 1.0e-4;
    let mut normal = [[0.0; 10]; 10];
    let mut rhs_x = [0.0; 10];
    let mut rhs_y = [0.0; 10];
    let mut rhs_z = [0.0; 10];
    for (row, &t) in params.iter().enumerate() {
        let basis = bernstein9(t);
        for i in 0..10 {
            rhs_x[i] += basis[i] * targets[row].x;
            rhs_y[i] += basis[i] * targets[row].y;
            rhs_z[i] += basis[i] * targets[row].z;
            for j in 0..10 {
                normal[i][j] += basis[i] * basis[j];
            }
        }
    }
    for i in 0..10 {
        normal[i][i] += lambda;
        rhs_x[i] += lambda * prior[i].x;
        rhs_y[i] += lambda * prior[i].y;
        rhs_z[i] += lambda * prior[i].z;
    }
    let xs = solve_10(normal, rhs_x);
    let ys = solve_10(normal, rhs_y);
    let zs = solve_10(normal, rhs_z);
    for index in 0..10 {
        set_boundary_control(
            patch,
            side,
            index,
            Point3::new(xs[index], ys[index], zs[index]),
        );
    }
}

fn boundary_controls(
    patch: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    side: BoundarySide,
) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = match side {
            BoundarySide::UMin => patch.grid[0][index],
            BoundarySide::UMax => patch.grid[9][index],
            BoundarySide::VMin => patch.grid[index][0],
            BoundarySide::VMax => patch.grid[index][9],
        };
    }
    out
}

fn trim_gap(
    patch_a: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    patch_b: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    trim: FullBoundaryTrim,
) -> f64 {
    let mut worst = 0.0_f64;
    for index in 0..9 {
        let fraction = index as f64 / 8.0;
        let param_a = lerp_scalar(trim.start_a, trim.end_a, fraction);
        let param_b = lerp_scalar(trim.start_b, trim.end_b, 1.0 - fraction);
        worst = worst.max(
            eval_boundary(patch_a, trim.side_a, param_a).distance_to(&eval_boundary(
                patch_b,
                trim.side_b,
                param_b,
            )),
        );
    }
    worst
}

fn max_boundary_displacement(
    before: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    after: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    side: BoundarySide,
) -> f64 {
    let mut worst = 0.0_f64;
    for index in 0..10 {
        let t = index as f64 / 9.0;
        worst =
            worst.max(eval_boundary(before, side, t).distance_to(&eval_boundary(after, side, t)));
    }
    worst
}

fn bernstein9(t: f64) -> [f64; 10] {
    let omt = 1.0 - t;
    let mut out = [0.0; 10];
    for (i, value) in out.iter_mut().enumerate() {
        *value = binomial9(i) * t.powi(i as i32) * omt.powi((9 - i) as i32);
    }
    out
}

fn binomial9(i: usize) -> f64 {
    [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0][i]
}

fn solve_10(mut matrix: [[f64; 10]; 10], mut rhs: [f64; 10]) -> [f64; 10] {
    for pivot in 0..10 {
        let mut best = pivot;
        for row in pivot + 1..10 {
            if matrix[row][pivot].abs() > matrix[best][pivot].abs() {
                best = row;
            }
        }
        matrix.swap(pivot, best);
        rhs.swap(pivot, best);
        let scale = matrix[pivot][pivot];
        if scale.abs() <= 1.0e-14 {
            continue;
        }
        for column in pivot..10 {
            matrix[pivot][column] /= scale;
        }
        rhs[pivot] /= scale;
        for row in 0..10 {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            for column in pivot..10 {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    rhs
}

fn set_boundary_control(
    patch: &mut popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    side: BoundarySide,
    index: usize,
    point: Point3,
) {
    match side {
        BoundarySide::UMin => patch.grid[0][index] = point,
        BoundarySide::UMax => patch.grid[9][index] = point,
        BoundarySide::VMin => patch.grid[index][0] = point,
        BoundarySide::VMax => patch.grid[index][9] = point,
    }
}

fn average_point(a: Point3, b: Point3) -> Point3 {
    a + (b - a) * 0.5
}

fn lerp_scalar(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

fn face_outer_half_edges(brep: &BRep, face: FaceHandle) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| {
            CoreError::InvalidHandle(
                "kp20 support-cell boundary reconciliation face invalid".into(),
            )
        })?
        .outer_loop;
    let start = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| {
            CoreError::InvalidHandle(
                "kp20 support-cell boundary reconciliation loop invalid".into(),
            )
        })?
        .half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Ok(out);
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 support-cell boundary reconciliation loop did not close".into(),
    ))
}
