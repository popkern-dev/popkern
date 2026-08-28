//! Shared test helpers for the apex_construction end-to-end test
//! suites.  Used by [`super::plan_e2e_tests`] and
//! [`super::plan_e2e_extras_tests`].
//!
//! Split out so each test file can stay under the 500-line
//! structural cap.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use crate::analysis;
use crate::fillet_surface::FilletSurfaceResult;
use crate::geometry_compute::compute_fillet_geometry;
use crate::types::{CapKind as TypesCapKind, ContinuityLevel, FilletRadius, FilletSpec};

pub(super) fn test_tolerance() -> Tolerance {
    Tolerance::default()
        .with_linear(1e-3)
        .with_angular(1e-3)
        .with_parametric(1e-6)
}

/// Collect the (apex, [edge; N]) tuple for the lowest-index
/// vertex of `brep`.  Picks all incident edges via half-edge
/// origin scan.
pub(super) fn pick_apex_and_incident_edges(brep: &BRep) -> (VertexHandle, Vec<EdgeHandle>) {
    let apex = brep.vertices.keys().next().expect("brep has vertices");
    let mut edge_set: std::collections::BTreeSet<EdgeHandle> = std::collections::BTreeSet::new();
    for he in brep.half_edges.values() {
        if he.origin == apex {
            edge_set.insert(he.edge);
        }
    }
    (apex, edge_set.into_iter().collect())
}

/// Build a list of [`FilletSurfaceResult`]s from a list of
/// (edge, radius) pairs, using the real
/// [`compute_fillet_geometry`] pipeline.
pub(super) fn build_fillet_surfaces(
    brep: &BRep,
    edges: &[EdgeHandle],
    radius: f64,
) -> Vec<FilletSurfaceResult> {
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let tol = test_tolerance();
    let analyses = analysis::analyze_edges(brep, edges, false)
        .expect("analyze_edges should succeed on uniform-radius corner");
    let mut surfaces = Vec::with_capacity(specs.len());
    for (spec, ana) in specs.iter().zip(analyses.iter()) {
        let geom = compute_fillet_geometry(
            ana,
            &spec.radius,
            brep,
            &tol,
            ContinuityLevel::G1,
            TypesCapKind::Open,
            crate::types::BlendProfile::Circular,
        )
        .expect("compute_fillet_geometry should succeed on cube/pyramid corner");
        surfaces.push(FilletSurfaceResult {
            surface: geom.surface.clone(),
            trim_curve_1: geom.trim_3d_1.clone(),
            trim_curve_2: geom.trim_3d_2.clone(),
            continuity_downgrade: None,
        });
    }
    surfaces
}

/// Like [`build_fillet_surfaces`] but each edge gets its own
/// radius (parallel slices over `edges` and `radii`).  Used by
/// variable-radius end-to-end tests.
pub(super) fn build_fillet_surfaces_per_edge(
    brep: &BRep,
    edges: &[EdgeHandle],
    radii: &[f64],
) -> Vec<FilletSurfaceResult> {
    assert_eq!(
        edges.len(),
        radii.len(),
        "edges and radii must have equal length"
    );
    let tol = test_tolerance();
    let analyses = analysis::analyze_edges(brep, edges, false)
        .expect("analyze_edges should succeed on variable-radius corner");
    let mut surfaces = Vec::with_capacity(edges.len());
    for ((edge, &r), ana) in edges.iter().zip(radii).zip(analyses.iter()) {
        let _ = edge;
        let radius_spec = FilletRadius::Constant(r);
        let geom = compute_fillet_geometry(
            ana,
            &radius_spec,
            brep,
            &tol,
            ContinuityLevel::G1,
            TypesCapKind::Open,
            crate::types::BlendProfile::Circular,
        )
        .expect("compute_fillet_geometry should succeed on variable-radius corner");
        surfaces.push(FilletSurfaceResult {
            surface: geom.surface.clone(),
            trim_curve_1: geom.trim_3d_1.clone(),
            trim_curve_2: geom.trim_3d_2.clone(),
            continuity_downgrade: None,
        });
    }
    surfaces
}

/// Replace `n` fillet surfaces (starting at index 0) with a
/// degenerate `Plane` variant that the planner's strict pass
/// rejects via `UnsupportedFilletKind` and the lenient pass
/// skips.  Trim curves are preserved so `compute_blend_boundary`
/// still succeeds.
pub(super) fn force_unsupported_kind(surfaces: &mut [FilletSurfaceResult], n: usize) {
    use popkern_core::geometry::Vec3;
    for fs in surfaces.iter_mut().take(n) {
        fs.surface = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
    }
}
