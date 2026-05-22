use super::boundary::PreparedExactBoundaryCurve;
use super::corner::{
    GeneratedInnerSpoke, PreparedExactBoundaryInnerSpokes, PreparedExactBoundarySectorBoundary,
};
use super::ribbon::{
    ExactBoundaryRibbonSample, PreparedExactBoundaryRibbon, PreparedExactBoundaryRibbons,
};
use super::sector::ExactBoundaryOuterSide;
use super::solve::{
    ExactBoundarySolveBudget, ExactBoundarySolveDeclineReason, solve_exact_boundary_surfaces,
};
use super::support::ExactBoundarySupportSurfaceFamily;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

#[test]
fn solve_locks_boundary_rows_and_derivative_controls() {
    let spokes = fixture_spokes(false);
    let ribbons = fixture_ribbons(spokes.sectors[0].fillet_boundary.edge);

    let solved =
        solve_exact_boundary_surfaces(&spokes, &ribbons, ExactBoundarySolveBudget::default())
            .unwrap();

    let sector = &spokes.sectors[0];
    let surface = &solved.surfaces[0].surface;
    assert_eq!(
        surface.num_u,
        sector.fillet_boundary.curve.control_points.len()
    );
    assert_eq!(
        surface.num_v,
        sector.face_boundary.curve.control_points.len()
    );
    for i in 1..surface.num_u - 1 {
        assert_eq!(
            surface.control_points[index(i, 0, surface.num_v)],
            sector.fillet_boundary.curve.control_points[i]
        );
        assert_eq!(
            surface.control_points[index(i, surface.num_v - 1, surface.num_v)],
            point_on_line(
                sector.face_inner_spoke.anchor,
                sector.face_inner_spoke.pole,
                i
            )
        );
    }
    for j in 0..surface.num_v {
        assert_eq!(
            surface.control_points[index(0, j, surface.num_v)],
            sector.face_boundary.curve.control_points[j]
        );
        assert_eq!(
            surface.control_points[index(surface.num_u - 1, j, surface.num_v)],
            point_on_line(
                sector.fillet_inner_spoke.anchor,
                sector.fillet_inner_spoke.pole,
                j
            )
        );
    }

    let derivative = Vec3::new(0.0, 0.2, 0.1);
    for i in 2..surface.num_u - 1 {
        let base = surface.control_points[index(i, 0, surface.num_v)];
        let next = surface.control_points[index(i, 1, surface.num_v)];
        assert!((next - base - derivative * (1.0 / 3.0)).length() <= 1.0e-12);
    }
    let face_derivative = Vec3::new(0.2, 0.0, 0.1);
    for j in 2..surface.num_v - 1 {
        let base = surface.control_points[index(0, j, surface.num_v)];
        let next = surface.control_points[index(1, j, surface.num_v)];
        assert!((next - base - face_derivative * (1.0 / 3.0)).length() <= 1.0e-12);
    }
    let fillet_twist =
        surface.control_points[index(1, 0, surface.num_v)] + derivative * (1.0 / 3.0);
    let face_twist =
        surface.control_points[index(0, 1, surface.num_v)] + face_derivative * (1.0 / 3.0);
    let averaged_twist = fillet_twist + (face_twist - fillet_twist) * 0.5;
    assert_eq!(
        surface.control_points[index(1, 1, surface.num_v)],
        averaged_twist
    );
}

#[test]
fn solve_is_deterministic_for_identical_inputs() {
    let spokes = fixture_spokes(false);
    let ribbons = fixture_ribbons(spokes.sectors[0].fillet_boundary.edge);

    let first =
        solve_exact_boundary_surfaces(&spokes, &ribbons, ExactBoundarySolveBudget::default())
            .unwrap();
    let second =
        solve_exact_boundary_surfaces(&spokes, &ribbons, ExactBoundarySolveBudget::default())
            .unwrap();

    assert_eq!(
        first.surfaces[0].surface.control_points,
        second.surfaces[0].surface.control_points
    );
    assert_eq!(
        first.surfaces[0].surface.knots_u,
        second.surfaces[0].surface.knots_u
    );
    assert_eq!(
        first.surfaces[0].surface.knots_v,
        second.surfaces[0].surface.knots_v
    );
}

#[test]
fn corner_conflict_declines_with_stable_diagnostic() {
    let spokes = fixture_spokes(true);
    let ribbons = fixture_ribbons(spokes.sectors[0].fillet_boundary.edge);

    let diagnostic =
        solve_exact_boundary_surfaces(&spokes, &ribbons, ExactBoundarySolveBudget::default())
            .expect_err("mismatched shared corner must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundarySolveDeclineReason::CornerDerivativeConflict
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.corner_derivative_conflict"
    );
}

#[test]
fn gross_foldover_gate_declines_collapsed_normals() {
    let spokes = fixture_spokes(false);
    let ribbons = fixture_ribbons(spokes.sectors[0].fillet_boundary.edge);

    let diagnostic = solve_exact_boundary_surfaces(
        &spokes,
        &ribbons,
        ExactBoundarySolveBudget {
            min_normal_length: 1.0e9,
            ..ExactBoundarySolveBudget::default()
        },
    )
    .expect_err("oversized normal gate must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundarySolveDeclineReason::GrossFoldover
    );
    assert_eq!(diagnostic.code, "exact_boundary_cap.solve.gross_foldover");
}

fn fixture_spokes(corner_conflict: bool) -> PreparedExactBoundaryInnerSpokes {
    let edge = dummy_edge();
    let shared = Point3::new(0.0, 0.0, 0.0);
    let fillet_anchor = Point3::new(1.0, 0.0, 0.0);
    let face_start = if corner_conflict {
        Point3::new(0.01, 0.0, 0.0)
    } else {
        shared
    };
    let face_anchor = Point3::new(0.0, 1.0, 0.0);
    let pole = Point3::new(1.0, 1.0, 0.2);
    let fillet_boundary = prepared_boundary(
        ExactBoundaryOuterSide::FilletEnd,
        edge,
        cubic_line(shared, fillet_anchor),
    );
    let face_boundary = prepared_boundary(
        ExactBoundaryOuterSide::FaceEdge,
        edge,
        cubic_line(face_start, face_anchor),
    );
    PreparedExactBoundaryInnerSpokes {
        sectors: vec![PreparedExactBoundarySectorBoundary {
            sector_index: 0,
            sector_name: "direct_apex.cap_sector.0".to_string(),
            fillet_boundary,
            face_boundary,
            fillet_inner_spoke: generated_spoke(
                "direct_apex.cap_sector.0.inner_spoke.fillet_anchor_to_pole",
                ExactBoundaryOuterSide::FilletEnd,
                fillet_anchor,
                pole,
            ),
            face_inner_spoke: generated_spoke(
                "direct_apex.cap_sector.0.inner_spoke.face_anchor_to_pole",
                ExactBoundaryOuterSide::FaceEdge,
                face_anchor,
                pole,
            ),
            outer_corner_gap: 0.0,
        }],
        pole,
    }
}

fn fixture_ribbons(edge: EdgeHandle) -> PreparedExactBoundaryRibbons {
    let support_half_edge = dummy_half_edge();
    let support_face = dummy_face();
    PreparedExactBoundaryRibbons {
        ribbons: vec![
            ribbon(
                ExactBoundaryOuterSide::FilletEnd,
                edge,
                support_half_edge,
                support_face,
                Vec3::new(0.0, 0.2, 0.1),
            ),
            ribbon(
                ExactBoundaryOuterSide::FaceEdge,
                edge,
                support_half_edge,
                support_face,
                Vec3::new(0.2, 0.0, 0.1),
            ),
        ],
    }
}

fn ribbon(
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    support_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
    cross_derivative: Vec3,
) -> PreparedExactBoundaryRibbon {
    PreparedExactBoundaryRibbon {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        side,
        edge,
        support_half_edge,
        support_face,
        support_family: ExactBoundarySupportSurfaceFamily::Plane,
        support_parameter_reversed: false,
        support_half_edge_curve_forward: true,
        max_support_residual: 0.0,
        samples: (0..4)
            .map(|i| {
                let parameter = i as f64 / 3.0;
                ExactBoundaryRibbonSample {
                    parameter,
                    support_uv: (parameter, 0.0),
                    boundary_point: Point3::new(parameter, 0.0, 0.0),
                    support_point: Point3::new(parameter, 0.0, 0.0),
                    boundary_tangent: Vec3::new(1.0, 0.0, 0.0),
                    support_normal: Vec3::new(0.0, 0.0, 1.0),
                    cross_derivative,
                }
            })
            .collect(),
    }
}

fn prepared_boundary(
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    curve: NurbsCurve,
) -> PreparedExactBoundaryCurve {
    PreparedExactBoundaryCurve {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        side,
        edge,
        curve,
    }
}

fn generated_spoke(
    name: &str,
    source_side: ExactBoundaryOuterSide,
    anchor: Point3,
    pole: Point3,
) -> GeneratedInnerSpoke {
    GeneratedInnerSpoke {
        name: name.to_string(),
        source_side,
        anchor,
        pole,
        curve: NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![anchor, pole],
            weights: vec![1.0, 1.0],
        },
    }
}

fn cubic_line(start: Point3, end: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: (0..4).map(|i| point_on_line(start, end, i)).collect(),
        weights: vec![1.0; 4],
    }
}

fn point_on_line(start: Point3, end: Point3, index: usize) -> Point3 {
    start + (end - start) * (index as f64 / 3.0)
}

fn index(i: usize, j: usize, num_v: usize) -> usize {
    i * num_v + j
}

fn dummy_edge() -> EdgeHandle {
    BRep::unit_cube(Tolerance::default())
        .unwrap()
        .edges
        .keys()
        .next()
        .unwrap()
}

fn dummy_half_edge() -> HalfEdgeHandle {
    BRep::unit_cube(Tolerance::default())
        .unwrap()
        .half_edges
        .keys()
        .next()
        .unwrap()
}

fn dummy_face() -> FaceHandle {
    BRep::unit_cube(Tolerance::default())
        .unwrap()
        .faces
        .keys()
        .next()
        .unwrap()
}
