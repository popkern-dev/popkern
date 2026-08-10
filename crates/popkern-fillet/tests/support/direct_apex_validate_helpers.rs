use popkern_core::error::ValidationError;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use popkern_fillet::error::BlendWarning;
use popkern_fillet::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{BlendResult, fillet_edges};

use crate::helpers::{incident_edges, pick_vertex_with_valence};

pub(crate) struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    pub(crate) fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        unsafe {
            std::env::set_var(key, value);
        }
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        unsafe {
            match &self.prior {
                Some(v) => std::env::set_var(self.key, v),
                None => std::env::remove_var(self.key),
            }
        }
    }
}

#[derive(Debug, Default)]
pub(crate) struct MetricBreakdown {
    pub(crate) structural: usize,
    pub(crate) path_a_relevant: usize,
    pub(crate) other_geometric: usize,
    pub(crate) worst_path_a_gap: f64,
}

pub(crate) fn breakdown(errors: &[ValidationError]) -> MetricBreakdown {
    let mut out = MetricBreakdown::default();
    for err in errors {
        match err {
            ValidationError::EulerPoincare { .. }
            | ValidationError::DanglingReference { .. }
            | ValidationError::BrokenTwin { .. }
            | ValidationError::BrokenCycle { .. }
            | ValidationError::InconsistentVertex { .. }
            | ValidationError::NonManifoldEdge { .. }
            | ValidationError::LowValence { .. } => out.structural += 1,
            ValidationError::OpenLoop { gap, .. } => {
                out.path_a_relevant += 1;
                if *gap > out.worst_path_a_gap {
                    out.worst_path_a_gap = *gap;
                }
            }
            ValidationError::ShellGap { gap, .. } => {
                out.path_a_relevant += 1;
                if *gap > out.worst_path_a_gap {
                    out.worst_path_a_gap = *gap;
                }
            }
            ValidationError::InconsistentOrientation { .. }
            | ValidationError::WindingDishonest { .. }
            | ValidationError::CurveSurfaceDeviation { .. }
            | ValidationError::TwistedSurface { .. }
            | ValidationError::Precision { .. }
            | ValidationError::SelfIntersection { .. } => out.other_geometric += 1,
        }
    }
    out
}

pub(crate) fn build_under_direct_apex(is_pyramid: bool, radii: &[f64]) -> BlendResult {
    build_under_direct_apex_with_pyramid_height(is_pyramid, radii, 0.8)
}

pub(crate) fn build_under_direct_apex_with_pyramid_height(
    is_pyramid: bool,
    radii: &[f64],
    pyramid_height: f64,
) -> BlendResult {
    let _g = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let tol = Tolerance::default();
    let brep = if is_pyramid {
        BRep::n_gonal_pyramid(radii.len(), 1.0, pyramid_height, tol).unwrap()
    } else {
        BRep::unit_cube(tol).unwrap()
    };
    let vertex = if is_pyramid {
        pick_vertex_with_valence(&brep, radii.len()).unwrap()
    } else {
        brep.vertices.keys().next().unwrap()
    };
    let edges: Vec<EdgeHandle> = incident_edges(&brep, vertex);
    assert_eq!(edges.len(), radii.len(), "edge count must match radii");
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        use_direct_apex_construction: true,
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        // These matrix fixtures pin native CONSTRUCTION machinery.
        // Their crossing-flagged serves (PK-245 triage debt) would be
        // re-routed to legacy by the honest serving-retry exclusion
        // (PK-258) — serving policy has its own pins; keep the retry
        // out of the construction measurement.
        serving_self_intersection_retry: false,
        ..FilletOptions::default()
    };
    fillet_edges(&brep, &specs, &options, tol, None)
        .expect("fillet_edges succeeds under DirectApex")
}

pub(crate) fn count_direct_apex_fallbacks(result: &BlendResult) -> usize {
    result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::DirectApexFellBackToLegacy { .. }))
        .count()
}

pub(crate) fn direct_apex_fallback_reasons(result: &BlendResult) -> Vec<String> {
    result
        .edge_warnings
        .iter()
        .filter_map(|(_, w)| match w {
            BlendWarning::DirectApexFellBackToLegacy { reason, .. } => Some(reason.clone()),
            _ => None,
        })
        .collect()
}

pub(crate) fn assert_native_gregory_cube_case(name: &str, radii: &[f64]) {
    assert_native_gregory_case(name, false, radii, 3);
}

/// PK-216 shared-corner containment pin: a UNIFORM repeated apex
/// (coincident sector corners) must decline the constructive Gregory
/// route with the class NAMED and serve the recorded legacy fallback
/// as a clean body — no glued surfaceless residuals, no dead sides.
/// History: the pre-PK-216 "native, 0 fallbacks" greens on these
/// fixtures were vacuous — the accepted partition body glued a
/// surfaceless duplicate apex octagon membrane over the strip ends
/// (mesh hole) plus a mis-scaled sector pillow, and none of the
/// gates censused it. Bisect-verified as old as the constructive
/// path (PK-125 already served it). Native quality returns with the
/// shared-corner substrate (PK-219).
/// Permanent postconditions on the welded native serve (PK-219):
/// the PK-216 membrane class (surfaceless residuals, dead boundary
/// sides) must never reach serving again.
pub(crate) fn assert_welded_native_postconditions(name: &str, result: &BlendResult) {
    let reasons = direct_apex_fallback_reasons(result);
    assert!(
        reasons.is_empty(),
        "{name} must serve natively through the welded shared-corner substrate; got {}:\n{}",
        reasons.len(),
        reasons.join("\n"),
    );
    let surfaceless = result
        .brep
        .faces
        .iter()
        .filter(|(_, f)| f.surface.is_none())
        .count();
    assert_eq!(
        surfaceless, 0,
        "{name} must not serve glued surfaceless residuals (the PK-216 membrane class)"
    );
    assert!(
        result.brep.boundary_edges.is_empty(),
        "{name} must serve without dead boundary sides"
    );
    let errors = result.brep.validate_all();
    let m = breakdown(&errors);
    assert_eq!(m.structural, 0, "{name} structural errors: {errors:#?}");
    assert_eq!(
        m.path_a_relevant, 0,
        "{name} open-loop/shell-gap errors: {errors:#?}"
    );
}

pub(crate) fn assert_native_gregory_case(
    name: &str,
    is_pyramid: bool,
    radii: &[f64],
    expected_distinct_cylinders: usize,
) {
    let result = build_under_direct_apex(is_pyramid, radii);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "{name} must complete through native direct-apex; got {} fallback(s):\n{}",
        reasons.len(),
        reasons
            .iter()
            .map(|r| format!("  - {r}"))
            .collect::<Vec<_>>()
            .join("\n"),
    );

    let errors = result.brep.validate_all();
    let m = breakdown(&errors);
    assert_eq!(
        m.structural, 0,
        "native {name} structural errors: {errors:#?}"
    );
    assert_eq!(
        m.path_a_relevant, 0,
        "native {name} open-loop/shell-gap errors: {errors:#?}"
    );

    assert_native_cap_has_cylinder_seams(name, &result, expected_distinct_cylinders);
}

pub(crate) fn assert_native_cap_has_cylinder_seams(
    name: &str,
    result: &BlendResult,
    expected_distinct_cylinders: usize,
) {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let cap_faces = native_cap_faces(result);
    let mut cylinder_seams = 0usize;
    let mut cylinder_faces = std::collections::BTreeSet::new();

    for (_, edge) in result.brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = result.brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = result.brep.half_edges.get(ha.twin) else {
            continue;
        };
        let Some(loop_a) = result.brep.loops.get(ha.loop_) else {
            continue;
        };
        let Some(loop_b) = result.brep.loops.get(hb.loop_) else {
            continue;
        };
        let face_a = loop_a.face;
        let face_b = loop_b.face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_faces.contains(&face_a);
        let b_is_cap = cap_faces.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let neighbor_face = if a_is_cap { face_b } else { face_a };
        let Some(neighbor_surface) = result
            .brep
            .faces
            .get(neighbor_face)
            .and_then(|f| f.surface.as_ref())
        else {
            continue;
        };
        if !matches!(neighbor_surface.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        cylinder_seams += 1;
        cylinder_faces.insert(neighbor_face);
    }

    assert!(
        cylinder_faces.len() >= expected_distinct_cylinders,
        "native {name} must expose at least {expected_distinct_cylinders} distinct total cap-to-cylinder \
         neighbor faces; got seams={cylinder_seams}, distinct_cylinders={}",
        cylinder_faces.len(),
    );
}

pub(crate) fn native_cap_faces(
    result: &BlendResult,
) -> std::collections::BTreeSet<popkern_core::topology::FaceHandle> {
    result
        .corner_patches
        .iter()
        .map(|p| p.face)
        .chain(
            result
                .repair_report
                .apex_routes
                .iter()
                .flat_map(|route| route.exact_boundary_report.committed_faces.iter().copied()),
        )
        .collect()
}
