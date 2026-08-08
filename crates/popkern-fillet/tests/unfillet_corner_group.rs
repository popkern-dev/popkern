//! PK-174 — recordless corner-group unfillet/modify.
//!
//! A cube corner filleted on its three edges serves three strips
//! meeting a corner cap. The record path (PK-143) replays the
//! operation; this matrix gates the RECORDLESS inverse: group Euler
//! surgery that restores the three sharp edges meeting at the
//! support-plane triple-intersection corner vertex.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn capped_corner_result(radii: [f64; 3]) -> popkern_fillet::BlendResult {
    capped_corner_result_with(radii, &FilletOptions::default())
}

fn capped_corner_result_with(
    radii: [f64; 3],
    options: &FilletOptions,
) -> popkern_fillet::BlendResult {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();
    let corner = brep
        .vertices
        .iter()
        .find(|(_, v)| {
            (v.position.x - 1.0).abs() < 1e-9
                && (v.position.y - 1.0).abs() < 1e-9
                && (v.position.z - 1.0).abs() < 1e-9
        })
        .map(|(vh, _)| vh)
        .expect("corner");
    let mut edges = Vec::new();
    for (e, d) in brep.edges.iter() {
        let he = &brep.half_edges[d.half_edge];
        let a = he.origin;
        let b = brep.half_edges[he.twin].origin;
        if a == corner || b == corner {
            edges.push(e);
        }
    }
    assert_eq!(edges.len(), 3);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii)
        .map(|(&edge, r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    fillet_edges(&brep, &specs, options, tol, None).expect("capped-corner fillet")
}

fn capped_corner(radii: [f64; 3]) -> (BRep, popkern_core::topology::FaceHandle) {
    let r = capped_corner_result(radii);
    // Pick a cylinder face that IDENTIFIES as a blend: the serving
    // body can also carry non-blend cylinder fragments (support
    // subdivision cells under the PK-167 in-place path), which the
    // identify layer correctly rejects.
    let strip = r
        .brep
        .faces
        .iter()
        .filter(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        })
        .map(|(fh, _)| fh)
        .find(|&fh| popkern_fillet::identify_blend(&r.brep, fh).is_some())
        .expect("identifiable strip face");
    (r.brep, strip)
}

/// Restored body must be exactly the unit cube: 6 faces, 8 vertices,
/// 12 edges, volume 1, validate_all clean, and the restored corner
/// vertex at (1,1,1).
fn assert_restored_cube(out: &BRep) {
    let errors = out.validate_all();
    assert!(errors.is_empty(), "restored cube must validate: {errors:?}");
    assert_eq!(
        (out.faces.len(), out.vertices.len(), out.edges.len()),
        (6, 8, 12),
        "restored cube census"
    );
    let corner_back = out.vertices.iter().any(|(_, v)| {
        (v.position.x - 1.0).abs() < 1e-6
            && (v.position.y - 1.0).abs() < 1e-6
            && (v.position.z - 1.0).abs() < 1e-6
    });
    assert!(
        corner_back,
        "triple-intersection corner vertex (1,1,1) must be restored"
    );
}

#[test]
fn recordless_uniform_corner_group_restores_cube() {
    let (brep, strip) = capped_corner([0.1, 0.1, 0.1]);
    eprintln!("[ucg] faces before: {}", brep.faces.len());
    for (fh, f) in brep.faces.iter() {
        let k = match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => "plane",
            Some(SurfaceGeometry::Cylinder { .. }) => "cylinder",
            Some(SurfaceGeometry::Sphere { .. }) => "sphere",
            Some(SurfaceGeometry::Nurbs(_)) => "nurbs",
            Some(SurfaceGeometry::Gregory(_)) => "gregory",
            Some(_) => "other",
            None => "none",
        };
        eprintln!("[ucg]   {fh:?}: {k}");
    }
    let out = popkern_fillet::unfillet_face(&brep, strip).expect("recordless corner-group restore");
    assert_restored_cube(&out);
}

/// The mixed-radius corner routes through the setback/Gregory
/// machinery and ships spanning sector faces plus far-corner
/// residuals — not a local corner group. The recordless restore
/// must DECLINE structurally (naming the record path) rather than
/// serve a partial body; `modify_fillet_radius_via_record` is the
/// commercial path for this class (PK-143).
#[test]
fn recordless_mixed_radius_setback_class_declines_structured() {
    let (brep, strip) = capped_corner([0.08, 0.1, 0.12]);
    let err = popkern_fillet::unfillet_face(&brep, strip)
        .expect_err("setback-class corner must decline recordless restore");
    let msg = format!("{err:?}");
    assert!(
        msg.contains("unfillet_via_record") || msg.contains("record"),
        "decline must route to the record path: {msg}"
    );
}

/// Recordless corner-group MODIFY: change one strip's radius on a
/// capped corner without a record. The group restores, then
/// re-fillets with the target at the new radius and the siblings at
/// their measured radii (the record path's coupled-radii semantics,
/// recordless).
#[test]
fn recordless_corner_group_modify_changes_target_radius() {
    let (brep, strip) = capped_corner([0.1, 0.1, 0.1]);
    let out = popkern_fillet::modify_fillet_radius(
        &brep,
        strip,
        popkern_fillet::FilletRadius::Constant(0.15),
        Tolerance::default(),
    )
    .expect("recordless corner-group modify");
    // The restored cube is exactly clean (surgery postcondition), so
    // residuals here are the FORWARD fillet's mixed-radius corner
    // debt: a direct (0.15,0.1,0.1) fillet of a pristine cube ships
    // the same single CurveSurfaceDeviation class (measured 1.4e-2;
    // see direct_mixed_15_10_10_fillet_baseline). Pin structural and
    // open cleanliness plus the single-defect ceiling; tighten to
    // empty when the mixed-corner class closes.
    let errors = out.validate_all();
    let non_deviation = errors
        .iter()
        .filter(|e| {
            !matches!(
                e,
                popkern_core::error::ValidationError::CurveSurfaceDeviation { .. }
            )
        })
        .count();
    assert_eq!(
        non_deviation, 0,
        "modified body must be structurally clean: {errors:?}"
    );
    assert!(
        errors.len() <= 1,
        "at most the known mixed-corner deviation: {errors:?}"
    );
    let mut radii: Vec<f64> = out
        .faces
        .iter()
        .filter_map(|(_, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Cylinder { radius, .. }) => Some(*radius),
            _ => None,
        })
        .collect();
    radii.sort_by(f64::total_cmp);
    assert_eq!(radii.len(), 3, "three strips after modify");
    assert!(
        (radii[0] - 0.1).abs() < 1e-9
            && (radii[1] - 0.1).abs() < 1e-9
            && (radii[2] - 0.15).abs() < 1e-9,
        "radii after modify: {radii:?}"
    );
}

/// Null-hypothesis probe for the modify gate: a direct fillet of a
/// pristine cube at the same mixed radii. Whatever deviation this
/// ships is the fillet machinery's mixed-corner class, not the
/// restore's.
#[test]
fn direct_mixed_15_10_10_fillet_baseline() {
    let (brep, _) = capped_corner([0.15, 0.1, 0.1]);
    let errors = brep.validate_all();
    eprintln!("[ucg-null] direct mixed fillet validate_all: {errors:?}");
    // Since the PK-245 default arming of the serving retry, the
    // DEFAULT serve for this corner re-routes to legacy (its
    // direct-apex candidate ships a strip×cap self-intersection)
    // and validates clean or, at worst, keeps the known
    // CurveSurfaceDeviation class. The retry contract itself is
    // pinned by `direct_mixed_15_10_10_retry_serves_clean`.
    assert!(
        errors.len() <= 1
            && errors.iter().all(|e| matches!(
                e,
                popkern_core::error::ValidationError::CurveSurfaceDeviation { .. }
            )),
        "mixed-corner default-serve debt grew beyond the pinned class: {errors:?}"
    );
}

/// PK-250 pin: with the serving-boundary retry opted in, the mixed
/// [0.15, 0.1, 0.1] corner's defective direct-apex candidate is
/// detected AFTER the full pipeline (strip 9v1 Cylinder r=0.1 × cap
/// sector Nurbs, contact 2.1e-7), the retry re-runs with the apex
/// excluded from direct-apex promotion, the legacy dispatch serves
/// the corner, and the served body validates clean — including the
/// self-intersection check. The retry is recorded on the repair
/// report with the offending face pair.
#[test]
fn direct_mixed_15_10_10_retry_serves_clean() {
    let options = FilletOptions {
        serving_self_intersection_retry: true,
        ..FilletOptions::default()
    };
    let r = capped_corner_result_with([0.15, 0.1, 0.1], &options);
    let retries = &r.repair_report.serving_self_intersection_retries;
    assert_eq!(
        retries.len(),
        1,
        "exactly one serving self-intersection retry: {retries:?}"
    );
    assert!(
        retries[0].detail.contains("Self-intersection"),
        "retry detail names the defect: {}",
        retries[0].detail
    );
    let excluded = retries[0].vertex;
    let route = r
        .repair_report
        .apex_routes
        .iter()
        .find(|route| route.vertex == excluded)
        .expect("excluded apex has a route record");
    assert!(
        !matches!(
            route.route,
            popkern_fillet::ApexConstructionRoute::DirectApex
        ),
        "excluded apex must not serve direct-apex on the retry: {:?}",
        route.route
    );
    let errors = r.brep.validate_all();
    assert!(
        errors.is_empty(),
        "retried serve must validate clean: {errors:?}"
    );
}
