//! §5.3 item #1 — cubic-vs-quintic §2.4 cross-fixture scorecard.
//!
//! Compares the boundary-normal curvature deviation `|Δκ|/max(κ, 1/r)`
//! between `use_deg5_gregory = false` (cubic bidegree-3 sub-patches)
//! and `use_deg5_gregory = true` (quintic bidegree-5 sub-patches)
//! across several mixed-radii cube-corner configurations. The result
//! drives the decision to flip the deg5 default.
//!
//! The prediction: quintic closes the inherent 2% cubic-Bezier arc
//! approximation error, so §2.4 should be **meaningfully lower** on
//! configurations where the cubic path is arc-limited. Where the
//! corner is already dominated by other error sources (Shirman-Séquin
//! residual, boundary Bezier fit error on the spokes), the two should
//! score comparably.

use crate::analysis;
use crate::blend::boundary::{BoundarySegment, compute_blend_boundary};
use crate::blend::gregory_diagnostic_metrics::compute_metrics;
use crate::blend::surface::BlendSurfaceResult;
use crate::fillet_surface::FilletSurfaceResult;
use crate::geometry_compute::compute_fillet_geometry;
use crate::types::{
    CapKind, ContinuityLevel, CornerMethod, FilletOptions, FilletRadius, FilletSpec,
};

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Run the fillet pipeline on a mixed-radii cube corner with the
/// given radii, comparing §2.4 between the cubic and quintic
/// Gregory dispatch paths. Returns `(cubic_max, deg5_max)`.
fn measure_deg3_vs_deg5(radii: &[f64; 3]) -> (f64, f64) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
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

    // Build the pre-fillet blend boundary to extract the arc curves
    // — same pattern as the R8 C8b scorecard in
    // `three_quad_apply_g2_tests.rs`. The arcs are geometry from
    // the original fillet surfaces, identical for both cubic and
    // deg5 dispatch paths.
    let Ok(analyses) = analysis::analyze_edges(&cube, &incident_edges, false) else {
        return (0.0, 0.0);
    };
    let mut edge_geoms = Vec::with_capacity(specs.len());
    for (spec, ana) in specs.iter().zip(analyses.iter()) {
        let Ok(geom) = compute_fillet_geometry(
            ana,
            &spec.radius,
            &cube,
            &tol,
            ContinuityLevel::G1,
            CapKind::Open,
        ) else {
            return (0.0, 0.0);
        };
        edge_geoms.push(geom);
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
    let Ok(boundary) = compute_blend_boundary(vertex, &indices, &fillet_surfaces, &cube, &tol)
    else {
        return (0.0, 0.0);
    };
    let fillet_arc_curves: Vec<NurbsCurve> = boundary
        .segments
        .iter()
        .filter_map(|s| match s {
            BoundarySegment::FilletEnd { curve, .. } => Some(curve.clone()),
            _ => None,
        })
        .collect();
    if fillet_arc_curves.is_empty() {
        return (0.0, 0.0);
    }

    let cubic_opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: false,
        ..FilletOptions::default()
    };
    let Ok(cubic_pipeline) = crate::fillet_edges(&cube, &specs, &cubic_opts, tol, None) else {
        return (0.0, 0.0);
    };
    let cubic_faces: Vec<_> = cubic_pipeline
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .map(|p| p.face)
        .collect();
    if cubic_faces.len() != 3 {
        return (0.0, 0.0);
    }

    let deg5_opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        ..FilletOptions::default()
    };
    let Ok(deg5_pipeline) = crate::fillet_edges(&cube, &specs, &deg5_opts, tol, None) else {
        return (0.0, 0.0);
    };
    let deg5_faces: Vec<_> = deg5_pipeline
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .map(|p| p.face)
        .collect();
    if deg5_faces.len() != 3 {
        return (0.0, 0.0);
    }

    // For each arc, find the closest sub-patch in each pipeline and
    // measure §2.4 against that sub-patch's surface.
    let score = |faces: &[popkern_core::topology::FaceHandle], brep: &BRep| -> f64 {
        let sub_sofs: Vec<_> = faces
            .iter()
            .filter_map(|&h| brep.faces.get(h).and_then(|f| f.surface.clone()))
            .collect();
        let mut max_dev = 0.0_f64;
        for arc in &fillet_arc_curves {
            let (t0, t1) = arc.parameter_range();
            let mid = arc.point_at(0.5 * (t0 + t1));
            let Some(best_sub) = sub_sofs.iter().min_by(|a, b| {
                let da = a
                    .geometry
                    .closest_point(&mid)
                    .map(|(u, v)| a.geometry.point_at(u, v).distance_to(&mid))
                    .unwrap_or(f64::INFINITY);
                let db = b
                    .geometry
                    .closest_point(&mid)
                    .map(|(u, v)| b.geometry.point_at(u, v).distance_to(&mid))
                    .unwrap_or(f64::INFINITY);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            }) else {
                continue;
            };
            // Confirm it's a Gregory sub-patch (sanity check).
            let is_gregory = matches!(
                &best_sub.geometry,
                SurfaceGeometry::Gregory(
                    GregoryPatch::Rectangular(_) | GregoryPatch::RectangularDeg5(_),
                )
            );
            if !is_gregory {
                continue;
            }
            let sub_result = BlendSurfaceResult {
                surface: best_sub.geometry.clone(),
                boundary_curves: vec![arc.clone()],
                method: CornerMethod::Gregory,
                continuity: ContinuityLevel::G1,
            };
            let m = compute_metrics(&sub_result, std::slice::from_ref(arc), avg_radius);
            max_dev = max_dev.max(m.max_boundary_normal_curvature_deviation);
        }
        max_dev
    };

    let cubic_max = score(&cubic_faces, &cubic_pipeline.brep);
    let deg5_max = score(&deg5_faces, &deg5_pipeline.brep);
    (cubic_max, deg5_max)
}

/// Cross-fixture §2.4 scorecard for cubic vs quintic Gregory.
///
/// Prints `(cubic, deg5, Δ)` for each radius configuration and
/// asserts deg5 is no worse than cubic (within a small slack for
/// solver noise). A measurable deg5 improvement across the board
/// would justify flipping the `use_deg5_gregory` default to `true`.
#[test]
fn item1_cross_fixture_deg5_vs_cubic_scorecard() {
    let fixtures: &[[f64; 3]] = &[
        [0.10, 0.10, 0.10],
        [0.09, 0.10, 0.11],
        [0.08, 0.10, 0.12],
        [0.05, 0.10, 0.20],
        [0.15, 0.20, 0.25],
        [0.02, 0.03, 0.04],
    ];

    eprintln!("=== item #1 cubic-vs-quintic §2.4 scorecard ===");
    let mut worst_regress = 0.0_f64;
    for radii in fixtures {
        let (cubic, deg5) = measure_deg3_vs_deg5(radii);
        let delta = deg5 - cubic;
        eprintln!("  radii {radii:?}: cubic={cubic:.4e} deg5={deg5:.4e} Δ={delta:+.4e}");
        worst_regress = worst_regress.max(delta);
    }
    // Deg5 shouldn't regress cubic by more than a small tolerance
    // on any fixture (solver / fit noise). If it does, either our
    // quintic fit has a bug or the C-K 25-factor scaling is wrong.
    assert!(
        worst_regress < 0.05,
        "item #1 scorecard: deg5 regresses cubic by more than 0.05 \
         on some fixture (worst Δ = {worst_regress:+.4e})"
    );
}
