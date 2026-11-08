//! Fixtures for the exact-boundary commit tests: solved surfaces and prepared pcurves.

use super::*;

pub(super) fn solved_surfaces() -> ExactBoundarySolvedSurfaces {
    ExactBoundarySolvedSurfaces {
        surfaces: vec![ExactBoundarySectorSurface {
            sector_index: 0,
            sector_name: "sector_0".to_string(),
            surface: unit_square_surface(),
            max_corner_gap: 0.0,
        }],
    }
}

pub(super) fn unit_square_surface() -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        weights: vec![1.0; 4],
        num_u: 2,
        num_v: 2,
    }
}

pub(super) fn prepared_pcurves() -> PreparedExactBoundaryPcurves {
    PreparedExactBoundaryPcurves {
        pairs: vec![
            pcurve_pair(ExactBoundaryOuterSide::FilletEnd),
            pcurve_pair(ExactBoundaryOuterSide::FaceEdge),
        ],
    }
}

pub(super) fn prepared_pcurves_for(seams: &SharedTrimSeams) -> PreparedExactBoundaryPcurves {
    let mut pcurves = prepared_pcurves();
    pcurves.pairs[0].edge = seams.edges[0];
    pcurves.pairs[0].support_half_edge = seams.support_half_edges[0];
    pcurves.pairs[0].support_pcurve = Pcurve::linear((1.0, 0.0), (0.0, 0.0), 0.0, 1.0);
    pcurves.pairs[1].edge = seams.edges[1];
    pcurves.pairs[1].support_half_edge = seams.support_half_edges[1];
    pcurves.pairs[1].cap_pcurve = Pcurve::linear((0.0, 0.25), (1.0, 0.25), 0.0, 1.0);
    pcurves.pairs[1].support_pcurve = Pcurve::linear((1.0, 0.25), (0.0, 0.25), 0.0, 1.0);
    pcurves
}

pub(super) fn pcurve_pair(side: ExactBoundaryOuterSide) -> PreparedExactBoundaryPcurvePair {
    PreparedExactBoundaryPcurvePair {
        sector_index: 0,
        sector_name: "sector_0".to_string(),
        side,
        edge: Default::default(),
        support_half_edge: HalfEdgeHandle::default(),
        support_face: FaceHandle::default(),
        support_surface: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        cap_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
        support_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
        max_cap_residual: 0.0,
        max_support_residual: 0.0,
        max_pair_residual: 0.0,
    }
}
