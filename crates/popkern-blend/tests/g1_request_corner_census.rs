//! PK-597: G1-request corner census gate. Every apex here once served a
//! creased or barely-G1 cap on a G1 request; each must now certify G1
//! over the whole blend (PK-499 certificates at the commercial gates),
//! served by the setback rung. The orthogonal cube corner is the
//! control: it keeps the exact rolling-ball sphere.

use std::collections::BTreeSet;

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{
    ApexConstructionDiagnostic, BlendResult, ContinuityLevel, CornerMethod, FilletOptions,
    FilletRadius, FilletSpec, SetbackFitMode, SetbackFitOutcome, fillet_edges,
};
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::healing::curvature::principal_curvatures;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

fn fillet_corner(
    body: &BRep,
    vertex: VertexHandle,
    radii: &[f64],
) -> (BlendResult, BTreeSet<FaceHandle>) {
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .into_iter()
        .zip(radii.iter().cycle())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G1,
        ..Default::default()
    };
    let pre = body.faces.keys().collect();
    let result = fillet_edges(body, &specs, &options, Tolerance::default(), None).expect("fillet");
    (result, pre)
}

/// Whole-blend G1 at the commercial gates, and the worst normal jump.
fn assert_certified_g1(
    name: &str,
    result: &BlendResult,
    pre: &BTreeSet<FaceHandle>,
    g1_bound: f64,
) {
    let blend: Vec<FaceHandle> = result
        .brep
        .faces
        .keys()
        .filter(|f| !pre.contains(f))
        .collect();
    let evidence = measure_blend_continuity_evidence(
        &result.brep,
        &blend,
        ContinuityLevel::G1,
        &ContinuityGates::commercial(),
        33,
    );
    assert_eq!(
        evidence.achieved,
        ContinuityLevel::G1,
        "{name}: {evidence:?}"
    );
    assert!(
        evidence.worst_g1_normal <= g1_bound,
        "{name}: worst g1 {:.3e} over {g1_bound:e}",
        evidence.worst_g1_normal
    );
    assert!(result.brep.validate().is_ok(), "{name}");
}

fn pyramid_apex(n: usize, h: f64, radii: &[f64]) -> (BlendResult, BTreeSet<FaceHandle>) {
    let body = BRep::n_gonal_pyramid(n, 1.0, h, Tolerance::default()).expect("pyramid");
    let apex = body
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex");
    fillet_corner(&body, apex, radii)
}

#[test]
fn oblique_uniform_apexes_certify_g1_on_the_setback_rung() {
    // Measured before PK-597: n = 4 served the boundary-ring Gregory cap
    // with a 27° cap-cap crease (body G0, g1 0.458). The rolling-ball
    // sphere the dispatch assumed owns the corner is never constructed at
    // an oblique apex; the setback rung now serves it (g1 9.3e-10).
    for (n, bound) in [(4, 1e-6), (5, 1e-6)] {
        let (result, pre) = pyramid_apex(n, 1.2, &[0.16]);
        assert!(
            result.corner_patches.iter().all(|p| p.setback_cap),
            "n={n}: setback rung serves: {:?}",
            result.corner_patches
        );
        assert_certified_g1(&format!("pyramid n={n} h=1.2 r=0.16"), &result, &pre, bound);
    }
}

#[test]
fn oblique_uniform_trihedral_apex_certifies_g1_on_the_setback_rung() {
    // PK-615: the oblique uniform n = 3 apex left the direct-apex route
    // (its exact-boundary cap pinches at the centre); the setback rung
    // serves it certified G1.
    let (result, pre) = pyramid_apex(3, 1.2, &[0.16]);
    assert!(
        result.corner_patches.iter().all(|p| p.setback_cap),
        "{:?}",
        result.corner_patches
    );
    assert_certified_g1("pyramid n=3 h=1.2 r=0.16", &result, &pre, 1e-6);
}

#[test]
fn mixed_radius_quad_apex_certifies_g1_on_the_setback_rung() {
    // Before PK-597 direct apex served the exact-boundary Coons cap with
    // cap-cap creases up to 28.6° (body G0, g1 0.479). Its served cap now
    // takes the crease re-run, and the setback rung serves: g1 2.8e-9.
    let (result, pre) = pyramid_apex(4, 0.8, &[0.04, 0.06, 0.05, 0.07]);
    assert!(
        result.corner_patches.iter().all(|p| p.setback_cap),
        "{:?}",
        result.corner_patches
    );
    assert_certified_g1("pyramid n=4 h=0.8 mixed", &result, &pre, 1e-6);
}

#[test]
fn orthogonal_cube_corner_keeps_the_rolling_ball_sphere() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(1.0, 1.0, 1.0)) < 1e-9)
        .map(|(h, _)| h)
        .expect("corner");
    let (result, pre) = fillet_corner(&cube, corner, &[0.2]);
    assert!(
        result
            .corner_patches
            .iter()
            .all(|p| p.method == CornerMethod::RollingBallSphere),
        "{:?}",
        result.corner_patches
    );
    assert_certified_g1("cube corner r=0.2", &result, &pre, 1e-6);
}

/// Largest sampled principal curvature over the corner caps' interiors
/// (30 x 30 per face), times `r_min`: the PK-612 fairness reading.
fn cap_fairness(result: &BlendResult, r_min: f64) -> f64 {
    let mut worst = 0.0_f64;
    for patch in &result.corner_patches {
        let Some(surface) = result
            .brep
            .faces
            .get(patch.face)
            .and_then(|f| f.surface.as_ref())
        else {
            continue;
        };
        let geometry: &dyn SurfaceEval = &surface.geometry;
        let ((u0, u1), (v0, v1)) = geometry.parameter_range();
        for i in 1..31 {
            for j in 1..31 {
                let u = u0 + (u1 - u0) * i as f64 / 31.0;
                let v = v0 + (v1 - v0) * j as f64 / 31.0;
                if let Some((k1, k2)) = principal_curvatures(geometry, u, v) {
                    worst = worst.max(k1.abs().max(k2.abs()));
                }
            }
        }
    }
    worst * r_min
}

/// A mixed-radius G1 row the corner-cell tangent fit rescues (PK-615):
/// served on the setback rung by `TangentCornerCells` after the plain
/// tangent fit failed, every corner piece claiming G1, the whole blend
/// certified G1 and the cap under the fairness bound.
fn assert_corner_cell_rescue(n: usize, h: f64, radii: &[f64]) {
    let name = format!("pyramid n={n} h={h} {radii:?}");
    let (result, pre) = pyramid_apex(n, h, radii);
    assert!(
        result
            .corner_patches
            .iter()
            .all(|p| p.setback_cap && p.achieved_continuity == ContinuityLevel::G1),
        "{name}: {:?}",
        result.corner_patches
    );
    let record = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|d| match d {
            ApexConstructionDiagnostic::SetbackFit(record) => Some(record),
            _ => None,
        })
        .unwrap_or_else(|| panic!("{name}: the rung's record"));
    let served = record
        .attempts
        .iter()
        .find(|a| matches!(a.outcome, SetbackFitOutcome::Served(_)))
        .unwrap_or_else(|| panic!("{name}: {record:#?}"));
    assert_eq!(served.mode, SetbackFitMode::TangentCornerCells, "{name}");
    assert_certified_g1(&name, &result, &pre, 1e-4);
    let r_min = radii.iter().take(n).copied().fold(f64::INFINITY, f64::min);
    let fairness = cap_fairness(&result, r_min);
    assert!(fairness < 4.0, "{name}: cap interior {fairness:.2} / r_min");
}

#[test]
fn mixed_radius_trihedral_g1_apexes_serve_the_corner_cell_tangent_fit() {
    // Measured 2026-10-01. Before PK-615 these served the boundary-ring
    // Gregory cap claiming G0 (pyramid-3 h0.4 wide: a 64° blend-blend
    // crease): the KP20 fan creased along a spoke (seam g1 0.79) and the
    // plain tangent fit met a rim at g1 0.99. With corner cells the
    // tangent fit certifies G1 at 2.1-3.2 / r_min.
    assert_corner_cell_rescue(3, 0.4, &[0.02, 0.04, 0.06]);
    assert_corner_cell_rescue(3, 0.4, &[0.05, 0.15, 0.08]);
    assert_corner_cell_rescue(3, 0.8, &[0.05, 0.15, 0.08]);
}

#[test]
fn mixed_radius_quad_g1_apexes_serve_the_corner_cell_tangent_fit() {
    // pyramid-4 h0.4 wide served the ring Gregory cap claiming G0 (fair
    // 102 / r_min); h0.8 variable the KP20 fan at 31.7 (PK-616's list).
    // Both now 1.0-1.4 / r_min.
    assert_corner_cell_rescue(4, 0.4, &[0.05, 0.15, 0.08, 0.2]);
    assert_corner_cell_rescue(4, 0.8, &[0.02, 0.04, 0.06, 0.04]);
}
