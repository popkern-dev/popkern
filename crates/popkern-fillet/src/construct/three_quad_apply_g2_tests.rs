//! §5.3 R8 C8b — direct §2.4 `|κ_n − 1/r|/max(κ_n, 1/r)` parity
//! measurement between polar-NURBS and three-quad Gregory paths
//! on a mixed-radius cube corner.
//!
//! Honest framing: both paths claim and deliver **G1** at the
//! blend boundary. §2.4's `< 0.01` is a **G2** threshold; G1
//! constructions aren't expected to hit it. What this test
//! guards is "rational-Gregory at G1 is no worse than polar-
//! NURBS at G1" — the R8 rewrite must not regress the G1
//! quality. True G2 (Chiyokura–Kimura twist + central-vertex
//! cross-patch solve) is C4b.2 work — now partially landed,
//! further refinement in C4b.2b.

use crate::analysis;
use crate::blend::boundary::{BoundarySegment, compute_blend_boundary};
use crate::blend::gregory_diagnostic_metrics::compute_metrics;
use crate::blend::surface::{BlendSurfaceResult, build_blend_surface};
use crate::fillet_surface::FilletSurfaceResult;
use crate::geometry_compute::compute_fillet_geometry;
use crate::types::{CapKind, ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, VertexMode};
use crate::vertex_classifier::classify_vertex_blend;

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn r8_c8b_boundary_normal_curvature_deviation_gregory_wins() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let radii = [0.08, 0.10, 0.12];
    let avg_radius = radii.iter().sum::<f64>() / 3.0;
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let tol = Tolerance::default();

    let analyses = analysis::analyze_edges(&cube, &incident_edges, false).unwrap();
    let mut edge_geoms = Vec::with_capacity(specs.len());
    for (spec, ana) in specs.iter().zip(analyses.iter()) {
        edge_geoms.push(
            compute_fillet_geometry(
                ana,
                &spec.radius,
                &cube,
                &tol,
                ContinuityLevel::G1,
                CapKind::Open,
            )
            .unwrap(),
        );
    }
    let fillet_surfaces: Vec<FilletSurfaceResult> = edge_geoms
        .iter()
        .map(|eg| FilletSurfaceResult {
            surface: eg.surface.clone(),
            trim_curve_1: eg.trim_3d_1.clone(),
            trim_curve_2: eg.trim_3d_2.clone(),
            continuity_downgrade: None,
        })
        .collect();
    let indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
    let boundary = compute_blend_boundary(vertex, &indices, &fillet_surfaces, &cube, &tol).unwrap();
    let fillet_arc_curves: Vec<NurbsCurve> = boundary
        .segments
        .iter()
        .filter_map(|s| match s {
            BoundarySegment::FilletEnd { curve, .. } => Some(curve.clone()),
            _ => None,
        })
        .collect();

    // Baseline: polar-NURBS single-surface corner blend.
    let classification = classify_vertex_blend(&cube, vertex, &incident_edges, &tol).classification;
    let polar_result = build_blend_surface(
        &boundary,
        &fillet_surfaces,
        &cube,
        VertexMode::Smooth,
        classification,
        &tol,
    )
    .expect("polar-NURBS blend surface must build");
    let polar_m = compute_metrics(&polar_result, &fillet_arc_curves, avg_radius);

    // Gregory (three-quad): run the full pipeline with the flag
    // on, then measure the §2.4 metric on each sub-patch.
    let gregory_opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    let gregory_pipeline =
        crate::fillet_edges(&cube, &specs, &gregory_opts, tol, None).expect("fillet_edges ok");
    let gregory_corner_faces: Vec<_> = gregory_pipeline
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .map(|p| p.face)
        .collect();
    assert_eq!(
        gregory_corner_faces.len(),
        3,
        "R8 C7 dispatch should yield 3 sub-patches at the corner"
    );

    let mut gregory_max_kappa_dev = 0.0_f64;
    for &face_h in &gregory_corner_faces {
        let face = gregory_pipeline.brep.faces.get(face_h).expect("sub-face");
        let sof = face.surface.as_ref().expect("sub-face surface");
        assert!(
            matches!(
                &sof.geometry,
                SurfaceGeometry::Gregory(GregoryPatch::Rectangular(_))
            ),
            "sub-face should hold a RectangularGregory surface"
        );
        let sub_result = BlendSurfaceResult {
            surface: sof.geometry.clone(),
            boundary_curves: fillet_arc_curves.clone(),
            method: crate::types::CornerMethod::Gregory,
            continuity: ContinuityLevel::G1,
        };
        let sub_m = compute_metrics(&sub_result, &fillet_arc_curves, avg_radius);
        gregory_max_kappa_dev =
            gregory_max_kappa_dev.max(sub_m.max_boundary_normal_curvature_deviation);
    }

    eprintln!(
        "R8 C8b §2.4 boundary-normal κ deviation scorecard (mixed radii 0.08/0.10/0.12):\n\
         polar-NURBS   max |Δκ|/max(κ,1/r) = {:.4e}  (target < 1e-2)\n\
         gregory (3Q)  max |Δκ|/max(κ,1/r) = {:.4e}  (target < 1e-2)\n\
         delta         = {:+.4e}",
        polar_m.max_boundary_normal_curvature_deviation,
        gregory_max_kappa_dev,
        gregory_max_kappa_dev - polar_m.max_boundary_normal_curvature_deviation,
    );

    assert!(
        gregory_max_kappa_dev <= polar_m.max_boundary_normal_curvature_deviation + 1e-6,
        "R8 C8b regression: Gregory boundary-κ deviation {gregory_max_kappa_dev:.4e} \
         exceeds polar-NURBS {:.4e}; expected ≤ polar for the flag flip to be safe",
        polar_m.max_boundary_normal_curvature_deviation,
    );
}
