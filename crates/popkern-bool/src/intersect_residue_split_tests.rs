//! Pins for the verdict-uniform residue split (PK-553): a residue that
//! runs through the other operand's wall is cut exactly on that wall;
//! a residue that stays on one side is left whole.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{OtherOperand, split_residues_at_verdict_changes};
use crate::intersect::{EndpointBoundaryType, IntersectionRecord, compute_face_aabbs};

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

/// A record on an A face of the cube (every record's `face_a`) along
/// `a -> b` in 3D, with owner-face uv `uv_a -> uv_b`; `residue` marks
/// an A-owned residue (judged against B, the same cube).
fn record(
    cube: &BRep,
    (a, b): (Point3, Point3),
    (uv_a, uv_b): (Point3, Point3),
    residue: bool,
) -> IntersectionRecord {
    let face = cube.faces.keys().next().expect("cube face");
    let uv = line(uv_a, uv_b);
    IntersectionRecord {
        face_a: face,
        face_b: face,
        curve: CurveOnSurface {
            curve_3d: line(a, b),
            curve_uv1: uv.clone(),
            curve_uv2: uv,
            is_closed: false,
        },
        is_tangent: false,
        is_coincident: false,
        start_boundary: EndpointBoundaryType::OnEdgeA,
        end_boundary: EndpointBoundaryType::OnEdgeA,
        residue_owner_a: residue.then_some(true),
    }
}

/// The cube's walls x = 0 and x = 1 as ordinary records on the owner
/// face, in its uv chart where u = x + 1.
fn walls(cube: &BRep) -> Vec<IntersectionRecord> {
    [0.0_f64, 1.0]
        .iter()
        .map(|&x| {
            record(
                cube,
                (Point3::new(x, 0.0, 0.37), Point3::new(x, 1.0, 0.37)),
                (
                    Point3::new(x + 1.0, -1.0, 0.0),
                    Point3::new(x + 1.0, 1.0, 0.0),
                ),
                false,
            )
        })
        .collect()
}

/// A residue along x in [x0, x1] at y = 0.5, z = 0.37 (u = x + 1).
fn residue_along(cube: &BRep, x0: f64, x1: f64) -> IntersectionRecord {
    record(
        cube,
        (Point3::new(x0, 0.5, 0.37), Point3::new(x1, 0.5, 0.37)),
        (
            Point3::new(x0 + 1.0, 0.0, 0.0),
            Point3::new(x1 + 1.0, 0.0, 0.0),
        ),
        true,
    )
}

fn split(cube: &BRep, residue: IntersectionRecord) -> Vec<IntersectionRecord> {
    let tol = Tolerance::default();
    let aabbs = compute_face_aabbs(cube);
    let mut records = walls(cube);
    records.push(residue);
    split_residues_at_verdict_changes(
        &mut records,
        &OtherOperand {
            brep_a: cube,
            brep_b: cube,
            aabbs_a: &aabbs,
            aabbs_b: &aabbs,
        },
        &tol,
    );
    records.retain(|r| r.residue_owner_a.is_some());
    records
}

#[test]
fn a_residue_through_the_other_operand_is_cut_on_its_walls() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    // In the cube's mid-plane, entering at x = 0 and leaving at x = 1.
    let pieces = split(&cube, residue_along(&cube, -1.0, 2.0));
    assert_eq!(pieces.len(), 3, "outside / inside / outside pieces");
    let ends: Vec<(f64, f64)> = pieces
        .iter()
        .map(|r| {
            let c = &r.curve.curve_3d;
            let (t0, t1) = c.parameter_range();
            (c.point_at(t0).x, c.point_at(t1).x)
        })
        .collect();
    for (got, want) in ends.iter().zip([(-1.0, 0.0), (0.0, 1.0), (1.0, 2.0)]) {
        assert!(
            (got.0 - want.0).abs() < 1e-9 && (got.1 - want.1).abs() < 1e-9,
            "piece {got:?} should span {want:?}"
        );
    }
    assert!(pieces.iter().all(|r| r.residue_owner_a == Some(true)));
}

#[test]
fn a_residue_on_one_side_of_the_other_operand_is_left_whole() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let outside = split(&cube, residue_along(&cube, 2.0, 3.0));
    assert_eq!(outside.len(), 1);
    let inside = split(&cube, residue_along(&cube, 0.2, 0.8));
    assert_eq!(inside.len(), 1);
}
