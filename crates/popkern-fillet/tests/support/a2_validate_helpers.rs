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
            | ValidationError::CurveSurfaceDeviation { .. }
            | ValidationError::TwistedSurface { .. }
            | ValidationError::Precision { .. }
            | ValidationError::SelfIntersection { .. } => out.other_geometric += 1,
        }
    }
    out
}

pub(crate) fn build_under_a2(is_pyramid: bool, radii: &[f64]) -> BlendResult {
    let _g = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let tol = Tolerance::default();
    let brep = if is_pyramid {
        BRep::n_gonal_pyramid(radii.len(), 1.0, 0.8, tol).unwrap()
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
        use_a2_apex_construction: true,
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    fillet_edges(&brep, &specs, &options, tol, None).expect("fillet_edges succeeds under A2")
}

pub(crate) fn count_a2_fallbacks(result: &BlendResult) -> usize {
    result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::A2FellBackToLegacy { .. }))
        .count()
}

pub(crate) fn a2_fallback_reasons(result: &BlendResult) -> Vec<String> {
    result
        .edge_warnings
        .iter()
        .filter_map(|(_, w)| match w {
            BlendWarning::A2FellBackToLegacy { reason, .. } => Some(reason.clone()),
            _ => None,
        })
        .collect()
}

pub(crate) fn assert_native_gregory_cube_case(name: &str, radii: &[f64]) {
    assert_native_gregory_case(name, false, radii, 3);
}

pub(crate) fn assert_native_gregory_case(
    name: &str,
    is_pyramid: bool,
    radii: &[f64],
    expected_distinct_cylinders: usize,
) {
    let result = build_under_a2(is_pyramid, radii);
    let reasons = a2_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "{name} must complete through native A2; got {} fallback(s):\n{}",
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

    let cap_faces = result
        .corner_patches
        .iter()
        .map(|p| p.face)
        .collect::<std::collections::BTreeSet<_>>();
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
        let face_a = result.brep.loops[ha.loop_].face;
        let face_b = result.brep.loops[hb.loop_].face;
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
