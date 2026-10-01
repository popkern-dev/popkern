//! KP cap fixtures for the STEP round-trip continuity audit: vertex
//! blends whose cap seams carry the kernel's hardest continuity
//! claims (the KP20 setback apex, a convex cube corner, and the
//! reflex saddle of a notched cube).

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::Point3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

/// The distinct edges incident to `vertex`, in handle order.
fn edges_at(body: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    edges
}

fn vertex_at(body: &BRep, at: Point3) -> Result<VertexHandle, String> {
    body.vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&at) < 1e-6)
        .map(|(h, _)| h)
        .ok_or_else(|| format!("no vertex at {at:?}"))
}

fn fillet_vertex(
    body: &BRep,
    vertex: VertexHandle,
    radii: &[f64],
    continuity: ContinuityLevel,
) -> Result<BRep, String> {
    let specs: Vec<FilletSpec> = edges_at(body, vertex)
        .into_iter()
        .zip(radii.iter().cycle())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity,
        ..Default::default()
    };
    fillet_edges(body, &specs, &options, Tolerance::default(), None)
        .map(|result| result.brep)
        .map_err(|e| format!("fillet_edges: {e:?}"))
}

/// The KP20 flagship: a pentagonal pyramid's five apex edges filleted
/// G2 at r = 0.16 (the setback cap serves the apex).
pub(crate) fn pyramid_apex() -> Result<BRep, String> {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).map_err(|e| format!("{e}"))?;
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .ok_or("no apex")?;
    fillet_vertex(&pyramid, apex, &[0.16], ContinuityLevel::G2)
}

/// A unit cube's (1,1,1) corner: three G2 edge fillets meeting in a
/// convex vertex blend.
pub(crate) fn cube_corner() -> Result<BRep, String> {
    let cube = BRep::unit_cube(Tolerance::default()).map_err(|e| format!("{e}"))?;
    let corner = vertex_at(&cube, Point3::new(1.0, 1.0, 1.0))?;
    fillet_vertex(&cube, corner, &[0.2], ContinuityLevel::G2)
}

/// A unit cube with its (+,+,+) corner notched out: the reflex saddle
/// vertex at (0.6, 0.6, 0.6) filleted with three radii (the saddle
/// setback cap).
pub(crate) fn notched_saddle() -> Result<BRep, String> {
    let tol = Tolerance::default();
    let notch = Transform3::translation(0.6, 0.6, 0.6).compose(&Transform3::scale_uniform(0.6));
    let cube = || BRep::unit_cube(tol).map_err(|e| format!("{e}"));
    let tool = cube()?.transform(&notch).map_err(|e| format!("{e}"))?;
    let body = popkern_bool::subtraction(cube()?, tool, tol, None).map_err(|e| format!("{e:?}"))?;
    let saddle = vertex_at(&body, Point3::new(0.6, 0.6, 0.6))?;
    fillet_vertex(&body, saddle, &[0.06, 0.08, 0.10], ContinuityLevel::G1)
}

/// The variational corner fit (PK-574) on the corner KP20 folds on: the
/// mixed-radius cube corner (0.08 / 0.2 / 0.3: the fan bends to
/// 468 / r_min), served by the curvature-weighted fit at certified G2
/// with no override. At 0.05 / 0.2 / 0.3 the G2 fit bends past the
/// fairness bound and the tangent fit serves (PK-614).
pub(crate) fn fit_cube_corner() -> Result<BRep, String> {
    let cube = BRep::unit_cube(Tolerance::default()).map_err(|e| format!("{e}"))?;
    let corner = vertex_at(&cube, Point3::new(1.0, 1.0, 1.0))?;
    fillet_vertex(&cube, corner, &[0.08, 0.2, 0.3], ContinuityLevel::G2)
}

/// The variational corner fit on the KP20 flagship with the KP20 cap
/// forced to decline: the tangent fit serves the pentagonal apex at
/// certified G1.
pub(crate) fn fit_pyramid_apex() -> Result<BRep, String> {
    popkern_blend::dispatch_overrides::with_kp20_cap_forced_to_decline(pyramid_apex)
}
