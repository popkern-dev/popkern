//! PK-174 — recordless corner-group unfillet/modify.
//!
//! A cube corner filleted on its three edges serves three strips
//! meeting a corner cap. The record path (PK-143) replays the
//! operation; this matrix gates the RECORDLESS inverse: group Euler
//! surgery that restores the three sharp edges meeting at the
//! support-plane triple-intersection corner vertex.

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn capped_corner_result(radii: [f64; 3]) -> popkern_blend::BlendResult {
    capped_corner_result_with(radii, &FilletOptions::default())
}

fn capped_corner_result_with(
    radii: [f64; 3],
    options: &FilletOptions,
) -> popkern_blend::BlendResult {
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
        .find(|&fh| popkern_blend::identify_blend(&r.brep, fh).is_some())
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
    let out = popkern_blend::unfillet_face(&brep, strip).expect("recordless corner-group restore");
    assert_restored_cube(&out);
}

/// Setback-served corners restore recordless (PK-608). The KP20 fan
/// (cube 0.08/0.1/0.12: 12 Gregory pieces) and the fitted cap (cube
/// 0.05/0.2/0.3: 3 NURBS sectors) leave more than one corner-to-far
/// spur after the group merge; before PK-608 the extra spurs skipped
/// two sector splits and the restore declined typed to the record path
/// (four edges 1.0 off a face folded over two support planes).
#[test]
fn recordless_setback_corner_groups_restore_cube() {
    for radii in [[0.08, 0.1, 0.12], [0.05, 0.2, 0.3]] {
        let result = capped_corner_result(radii);
        assert!(
            result.corner_patches.iter().all(|p| p.setback_cap),
            "{radii:?}: the setback rung owns the corner: {:?}",
            result.corner_patches
        );
        let (brep, strip) = capped_corner(radii);
        let out = popkern_blend::unfillet_face(&brep, strip)
            .unwrap_or_else(|e| panic!("{radii:?}: recordless restore: {e:?}"));
        assert_restored_cube(&out);
    }
}

/// Recordless MODIFY of a setback-served corner (PK-608): the group
/// restores, then re-fillets with the target at the new radius.
#[test]
fn recordless_setback_corner_group_modify_changes_target_radius() {
    let (brep, strip) = capped_corner([0.08, 0.1, 0.12]);
    let out = popkern_blend::modify_fillet_radius(
        &brep,
        strip,
        popkern_blend::FilletRadius::Constant(0.14),
        Tolerance::default(),
    )
    .expect("recordless setback corner modify");
    assert!(out.validate().is_ok(), "{:?}", out.validate());
    let radii: std::collections::BTreeSet<u64> = out
        .faces
        .iter()
        .filter_map(|(_, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Cylinder { radius, .. }) => Some((radius * 1e6).round() as u64),
            _ => None,
        })
        .collect();
    assert!(
        radii.contains(&140_000),
        "the target strip is re-filleted at 0.14: {radii:?}"
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
    let out = popkern_blend::modify_fillet_radius(
        &brep,
        strip,
        popkern_blend::FilletRadius::Constant(0.15),
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
    // DEFAULT serve for this corner re-routes to vertex-blend (its
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

/// PK-250/PK-253/PK-597 pin: with the serving-boundary retry opted in,
/// the mixed [0.15, 0.1, 0.1] corner's direct-apex candidate crosses a
/// strip (strip 9v1 Cylinder r=0.1 × cap sector Nurbs) and the PK-251
/// footprint repair clears the crossing in place — but the repaired cap
/// still creases across its own seams (certified cap-cap normal sine
/// 2.4e-1, body G0), so PK-597's crease check sends the corner to the
/// exclusion re-run. The retry record names the crease; the served body
/// validates clean and certifies G1 (g1 1.6e-6).
#[test]
fn direct_mixed_15_10_10_retry_serves_clean() {
    let options = FilletOptions {
        serving_self_intersection_retry: true,
        ..FilletOptions::default()
    };
    let r = capped_corner_result_with([0.15, 0.1, 0.1], &options);
    let repairs = &r.repair_report.serving_footprint_repairs;
    assert_eq!(
        repairs.len(),
        1,
        "exactly one serving footprint repair: {repairs:?}"
    );
    assert!(
        repairs[0].cleared,
        "the repair clears the crossing in place: {repairs:?}"
    );
    let retries = &r.repair_report.serving_self_intersection_retries;
    assert!(
        retries
            .iter()
            .any(|retry| retry.detail.contains("served cap-cap crease")),
        "the creased cap takes the exclusion re-run: {retries:?}"
    );
    assert!(
        r.corner_patches.iter().all(|p| p.setback_cap),
        "the setback rung serves the corner: {:?}",
        r.corner_patches
    );
    let errors = r.brep.validate_all();
    assert!(
        errors.is_empty(),
        "the re-run serve must validate clean: {errors:?}"
    );
}
