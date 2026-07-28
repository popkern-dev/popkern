//! Shared low-valence apex harness for the PnS gates and probes.
//!
//! Fillets an n-gonal pyramid apex and reports what the pipeline
//! actually committed: continuity metrics over the live seam
//! topology, the construction route that served the corner (read from
//! `CornerPatchInfo::method`, never re-derived), and the structured
//! `pns.*` diagnostics recorded along the way.
//!
//! Consumed via `#[path = "support/pns_apex_outcome.rs"]` by
//! `pns_low_valence_gates.rs` (the gates) and
//! `pns_low_valence_probes.rs` (the ignored diagnosis probes), so both
//! measure through the identical harness. The including crate must
//! also declare `mod kp_gate_helpers;` at its root.
//!
//! Each including binary reads a subset of `ApexOutcome`'s fields;
//! the rest are dead code from that binary's perspective, hence the
//! top-level allow (same idiom as `kp_gate_helpers.rs`).
#![allow(dead_code)]

use crate::kp_gate_helpers as helpers;
use helpers::{N_SAMPLES, incident_edges, measure_live_cap_seams, pick_vertex_with_valence};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::types::{CornerMethod, FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{ContinuityLevel, fillet_edges};

/// What a low-valence apex fillet produced: the continuity metrics
/// plus the construction route that actually served the corner.
pub struct ApexOutcome {
    pub metrics: helpers::LiveSeamMetrics,
    /// Distinct `CornerMethod`s across the committed corner patches.
    pub methods: Vec<CornerMethod>,
    /// Distinct achieved continuity levels across those patches.
    pub achieved: Vec<ContinuityLevel>,
    pub cap_faces: usize,
    /// The committed cap faces themselves, for anatomy probes.
    pub cap_face_handles: Vec<popkern_core::topology::FaceHandle>,
    /// The BRep the corner shipped in.
    pub brep: BRep,
    pub valid: bool,
    /// Structured `pns.*` diagnostics the pipeline recorded for this
    /// corner — the record of why PnS did or did not serve it.
    pub pns_diagnostics: Vec<String>,
    /// Per-seam breakdown: which surface kinds meet, and how badly.
    pub seam_detail: Vec<String>,
}

impl ApexOutcome {
    pub fn served_by_pns(&self) -> bool {
        self.methods.contains(&CornerMethod::PolyhedralNetSpline)
    }
}

/// Fillet an n-gonal pyramid apex and measure the committed cap
/// across its live seam topology, recording which route served it.
///
/// Route identity is read from the pipeline's own
/// `CornerPatchInfo::method`, not re-derived here — a gate that
/// guessed the route could pass while measuring a cap the route under
/// test never built.
pub fn run_pyramid_apex_outcome(n: usize, radii: &[f64]) -> ApexOutcome {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, 1.0, 3.0, Tolerance::default())
        .expect("n_gonal_pyramid fixture must build");
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    run_corner_outcome(pyr, vertex, radii)
}

/// Fillet one corner of the unit cube — PK-192's named `n = 3`
/// fixture class (three planes meeting at a trihedral corner) — and
/// measure through the identical instrument as the pyramid harness,
/// so the two fixtures' numbers are directly comparable.
pub fn run_cube_corner_outcome(radii: &[f64]) -> ApexOutcome {
    assert_eq!(radii.len(), 3, "a cube corner has exactly 3 edges");
    let cube = BRep::unit_cube(Tolerance::default()).expect("unit cube fixture must build");
    let vertex = cube.vertices.keys().next().expect("unit cube has vertices");
    run_corner_outcome(cube, vertex, radii)
}

/// Same pyramid harness with the KP-polar ladder disabled, so the
/// PnS dispatch inside `insert_vertex_blend` actually arbitrates at
/// `n >= 5` (with the ladder on, the setback/KP-polar route claims
/// those corners upstream and PnS never bids). Probe-only control for
/// valence-scaling comparisons.
pub fn run_pyramid_apex_outcome_no_polar(n: usize, radii: &[f64]) -> ApexOutcome {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, 1.0, 3.0, Tolerance::default())
        .expect("n_gonal_pyramid fixture must build");
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    run_corner_outcome_with(pyr, vertex, radii, false)
}

fn run_corner_outcome(
    brep: BRep,
    vertex: popkern_core::topology::VertexHandle,
    radii: &[f64],
) -> ApexOutcome {
    run_corner_outcome_with(brep, vertex, radii, true)
}

fn run_corner_outcome_with(
    brep: BRep,
    vertex: popkern_core::topology::VertexHandle,
    radii: &[f64],
    kp_polar: bool,
) -> ApexOutcome {
    let edges = incident_edges(&brep, vertex);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: kp_polar,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&brep, &specs, &opts, Tolerance::default(), None)
        .expect("low-valence corner fillet must not error");
    let mut methods: Vec<CornerMethod> = Vec::new();
    let mut achieved: Vec<ContinuityLevel> = Vec::new();
    for patch in &result.corner_patches {
        if !methods.contains(&patch.method) {
            methods.push(patch.method);
        }
        if !achieved.contains(&patch.achieved_continuity) {
            achieved.push(patch.achieved_continuity);
        }
    }
    let faces = helpers::cap_outer_faces_from_result(&result);
    let valid = result.brep.validate().is_ok();
    let metrics = measure_live_cap_seams(&result.brep, &faces, N_SAMPLES);
    let pns_diagnostics: Vec<String> = result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(e)
                if e.stage.contains("pns") =>
            {
                Some(format!(
                    "{}[{}]{}",
                    e.stage,
                    e.code,
                    e.detail_code
                        .as_ref()
                        .map(|d| format!(" {d}"))
                        .unwrap_or_default()
                ))
            }
            _ => None,
        })
        .collect();
    let seam_detail = describe_cap_seams(&result.brep, &faces);
    ApexOutcome {
        metrics,
        methods,
        achieved,
        cap_faces: faces.len(),
        cap_face_handles: faces.clone(),
        brep: result.brep.clone(),
        valid,
        pns_diagnostics,
        seam_detail,
    }
}

/// Name each measured cap-cap seam by the surface kinds it joins, so
/// a bad number can be attributed to a real pair of faces rather than
/// blamed on the measurement.
fn describe_cap_seams(
    brep: &popkern_core::topology::BRep,
    faces: &[popkern_core::topology::FaceHandle],
) -> Vec<String> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let kind = |f: popkern_core::topology::FaceHandle| -> String {
        match brep.faces.get(f).and_then(|d| d.surface.as_ref()) {
            Some(s) => match &s.geometry {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cylinder",
                SurfaceGeometry::Cone { .. } => "Cone",
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Torus { .. } => "Torus",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                SurfaceGeometry::Gregory(_) => "Gregory",
                SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
            }
            .to_string(),
            None => "none".to_string(),
        }
    };
    let set: std::collections::BTreeSet<_> = faces.iter().copied().collect();
    let mut out = Vec::new();
    for (_eh, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(twin.loop_)) else {
            continue;
        };
        let (fa, fb) = (la.face, lb.face);
        if fa == fb || !set.contains(&fa) || !set.contains(&fb) {
            continue;
        }
        let m = measure_live_cap_seams(brep, &[fa, fb], N_SAMPLES);
        let profile =
            helpers::seam_cos_profile(brep, fa, fb, &[0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95])
                .map(|p| {
                    p.iter()
                        .map(|(f, c)| format!("{f:.2}:{c:.4}"))
                        .collect::<Vec<_>>()
                        .join(" ")
                })
                .unwrap_or_else(|| "<no profile>".to_string());
        out.push(format!(
            "{} x {} -> g0={:.2e} g1={:.2e} cos_min={:.5}  cos(f)= {profile}",
            kind(fa),
            kind(fb),
            m.worst_g0_gap,
            m.worst_g1_normal,
            m.cos_min
        ));
    }
    out
}
