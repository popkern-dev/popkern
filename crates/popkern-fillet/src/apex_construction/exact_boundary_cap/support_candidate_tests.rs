use super::facade::{ExactBoundaryCapCandidateStatus, run_exact_boundary_cap_candidate};
use super::sector::{ExactBoundaryOuterSeam, extract_exact_boundary_cap_sectors};
use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{
    TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology,
};
use crate::apex_construction::telemetry::PlanTelemetry;
use crate::types::ApexConstructionDiagnostic;
use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, SurfaceEval, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};
use std::collections::BTreeMap;
use std::f64::consts::TAU;

#[test]
fn nurbs_supports_reach_full_candidate_with_family_performance() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors = trihedral_sectors(&edges);
    install_nurbs_plane_supports(&mut brep, cap_face, &plan, &sectors);

    let result = run_exact_boundary_cap_candidate(&mut brep, cap_face, &plan, &sectors)
        .expect("nurbs support candidate should produce a structured result");

    assert_eq!(result.performance.nurbs_support_ribbon_count, 6);
    assert_eq!(result.performance.support_ribbon_count, 6);
    assert!(result.performance.support_sample_count >= 30);
    assert!(result.performance.elapsed_micros > 0);
    assert!(
        !has_exact_boundary_code(
            &result.diagnostics,
            "exact_boundary_cap.ribbon.support_pcurve_stale"
        ),
        "trusted NURBS support fixture should not fail support pcurve readiness"
    );
}

#[test]
fn sphere_and_cone_full_candidates_decline_with_structured_support_diagnostics() {
    for surface in [
        SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        },
        SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: 0.4,
        },
    ] {
        let (mut brep, cap_face, edges) = hex_cap_fixture();
        let plan = trihedral_plan(&brep, &edges);
        let sectors = trihedral_sectors(&edges);
        let exact_sectors =
            extract_exact_boundary_cap_sectors(&brep, cap_face, &plan, &sectors).unwrap();
        let face = exact_sectors[0].outer_seams[0].support_face;
        brep.faces[face].surface = Some(SurfaceOnFace {
            geometry: surface.clone(),
            orientation: true,
        });

        let result = run_exact_boundary_cap_candidate(&mut brep, cap_face, &plan, &sectors)
            .expect("analytic support mismatch should produce a structured result");

        assert_eq!(result.status, ExactBoundaryCapCandidateStatus::Declined);
        assert!(result.performance.elapsed_micros > 0);
        assert!(
            has_exact_boundary_code(
                &result.diagnostics,
                "exact_boundary_cap.ribbon.support_pcurve_stale",
            ),
            "sphere/cone support mismatch must decline with a stable pcurve readiness diagnostic"
        );
    }
}

#[test]
fn periodic_cylinder_supports_reach_full_candidate_with_family_performance() {
    let (mut brep, cap_face, edges) = hex_cap_fixture();
    let plan = trihedral_plan(&brep, &edges);
    let sectors = trihedral_sectors(&edges);
    install_periodic_cylinder_supports(&mut brep, cap_face, &plan, &sectors);

    let result = run_exact_boundary_cap_candidate(&mut brep, cap_face, &plan, &sectors)
        .expect("periodic cylinder support candidate should produce a structured result");

    assert_eq!(result.performance.cylinder_support_ribbon_count, 6);
    assert_eq!(result.performance.support_ribbon_count, 6);
    assert!(result.performance.elapsed_micros > 0);
    assert!(
        !has_exact_boundary_code(
            &result.diagnostics,
            "exact_boundary_cap.ribbon.support_pcurve_stale"
        ),
        "periodic cylinder pcurve wrapping should not be classified as stale"
    );
}

fn install_nurbs_plane_supports(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    sectors: &[GregoryCapSector],
) {
    let exact_sectors = extract_exact_boundary_cap_sectors(brep, cap_face, plan, sectors).unwrap();
    let mut uv_ranges: BTreeMap<FaceHandle, ((f64, f64), (f64, f64))> = BTreeMap::new();
    for seam in exact_sectors
        .iter()
        .flat_map(|sector| sector.outer_seams.iter())
    {
        let pcurve = brep.half_edges[seam.support_half_edge]
            .pcurve
            .as_ref()
            .expect("fixture support pcurve exists");
        let start = pcurve.point_at(pcurve.t_start);
        let end = pcurve.point_at(pcurve.t_end);
        let entry = uv_ranges
            .entry(seam.support_face)
            .or_insert(((start.x, start.x), (start.y, start.y)));
        entry.0.0 = entry.0.0.min(start.x).min(end.x);
        entry.0.1 = entry.0.1.max(start.x).max(end.x);
        entry.1.0 = entry.1.0.min(start.y).min(end.y);
        entry.1.1 = entry.1.1.max(start.y).max(end.y);
    }
    for (face, ((u0, u1), (v0, v1))) in uv_ranges {
        let surface = brep.faces[face].surface.as_ref().unwrap().geometry.clone();
        let (u0, u1) = padded_range(u0, u1);
        let (v0, v1) = padded_range(v0, v1);
        brep.faces[face].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(bilinear_surface_from(&surface, u0, u1, v0, v1)),
            orientation: true,
        });
    }
}

fn install_periodic_cylinder_supports(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    sectors: &[GregoryCapSector],
) {
    let exact_sectors = extract_exact_boundary_cap_sectors(brep, cap_face, plan, sectors).unwrap();
    for seam in exact_sectors
        .iter()
        .flat_map(|sector| sector.outer_seams.iter())
    {
        install_periodic_cylinder_support(brep, seam);
    }
}

fn install_periodic_cylinder_support(brep: &mut BRep, seam: &ExactBoundaryOuterSeam) {
    let curve = brep.edges[seam.edge].curve.as_ref().unwrap();
    let start = curve.geometry.point_at(curve.t_start);
    let end = curve.geometry.point_at(curve.t_end);
    let axis_vector = end - start;
    let length = axis_vector.length();
    let axis = axis_vector.normalized().unwrap();
    let (radial, _) = orthonormal_basis(axis);
    brep.faces[seam.support_face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Cylinder {
            origin: start - radial,
            axis,
            radius: 1.0,
        },
        orientation: true,
    });
    brep.half_edges[seam.support_half_edge].pcurve =
        Some(Pcurve::linear((TAU, 0.0), (TAU, length), 0.0, 1.0));
}

fn bilinear_surface_from(
    surface: &SurfaceGeometry,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
) -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![u0, u0, u1, u1],
        knots_v: vec![v0, v0, v1, v1],
        control_points: vec![
            surface.point_at(u0, v0),
            surface.point_at(u0, v1),
            surface.point_at(u1, v0),
            surface.point_at(u1, v1),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    }
}

fn padded_range(a: f64, b: f64) -> (f64, f64) {
    if (b - a).abs() <= 1.0e-9 {
        (a - 1.0, b + 1.0)
    } else {
        (a, b)
    }
}

fn has_exact_boundary_code(diagnostics: &[ApexConstructionDiagnostic], code: &str) -> bool {
    diagnostics.iter().any(|diagnostic| {
        matches!(
            diagnostic,
            ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic) if diagnostic.code == code
        )
    })
}

fn hex_cap_fixture() -> (BRep, FaceHandle, Vec<EdgeHandle>) {
    let brep = BRep::n_gonal_pyramid(6, 1.0, 0.8, Tolerance::default()).unwrap();
    let cap_face = brep
        .faces
        .keys()
        .find(|&face| loop_edges(&brep, face).len() == 6)
        .unwrap();
    let edges = loop_edges(&brep, cap_face);
    (brep, cap_face, edges)
}

fn loop_edges(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut cursor = start;
    let mut edges = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = &brep.half_edges[cursor];
        edges.push(half_edge.edge);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    edges
}

fn trihedral_sectors(edges: &[EdgeHandle]) -> Vec<GregoryCapSector> {
    (0..3)
        .map(|sector_index| GregoryCapSector {
            sector_index,
            fillet_index: sector_index,
            fillet_end_edge: edges[2 * sector_index],
            face_edge: edges[2 * sector_index + 1],
        })
        .collect()
}

fn trihedral_plan(brep: &BRep, edges: &[EdgeHandle]) -> ApexPlan {
    ApexPlan {
        apex_vertex: brep.vertices.keys().next().unwrap(),
        cap_kind: CapKind::Gregory {
            per_fillet_terminations: vec![
                [Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)];
                3
            ],
            per_fillet_caps: vec![(Point3::new(0.0, 0.0, 0.0), 1.0); 3],
        },
        target_topology: TargetTopology {
            apex_vertices: Vec::new(),
            apex_edges: (0..6)
                .map(|index| TargetEdge {
                    origin_vertex_idx: index,
                    dest_vertex_idx: (index + 1) % 6,
                    curve: dummy_curve(),
                    role: if index % 2 == 0 {
                        TargetEdgeRole::CapVsCyl {
                            fillet_idx: index / 2,
                        }
                    } else {
                        TargetEdgeRole::CapVsPlane {
                            shared_face_idx: index / 2,
                        }
                    },
                })
                .collect(),
            cap_face: TargetCapFace::Gregory {
                outer_loop_edge_indices: (0..edges.len()).collect(),
            },
            face_boundary_updates: Vec::new(),
        },
        telemetry: PlanTelemetry::default(),
    }
}

fn dummy_curve() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }
}
