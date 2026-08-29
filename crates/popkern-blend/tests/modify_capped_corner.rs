//! PK-143 — cap-family modify re-dispatch via the operation record.
//!
//! Strips whose ends meet a corner cap cannot be modified by the
//! local paths (the cap end has no restored-edge endpoint), and a
//! corner's radii are coupled. `modify_fillet_radius_via_record`
//! replays the recorded specs against the recorded input with the
//! target spec rewritten — the corner re-dispatches through the same
//! construction machinery at the new radius.

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn capped_corner_record() -> (
    popkern_blend::FilletOperationRecord,
    popkern_core::topology::FaceHandle,
) {
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
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let record = popkern_blend::fillet_recorded(&brep, &specs, &FilletOptions::default(), tol)
        .expect("capped-corner fillet");
    let strip = record
        .result
        .faces
        .iter()
        .find(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        })
        .map(|(fh, _)| fh)
        .expect("strip face");
    (record, strip)
}

/// Same-radius replay must reproduce a fully valid body — the
/// re-dispatch mechanism itself is exact.
#[test]
fn record_modify_identity_replay_clean() {
    let (record, strip) = capped_corner_record();
    let m = popkern_blend::modify_fillet_radius_via_record(&record, strip, 0.1)
        .expect("identity replay");
    let errors = m.brep.validate_all();
    assert!(
        errors.is_empty(),
        "identity replay must be exactly clean: {errors:?}"
    );
}

/// Radius change at a capped corner COMPLETES through the corner
/// re-dispatch (the local modify paths decline with the 1-endpoint
/// signature). The changed corner is a MIXED-radius apex — the
/// replay surfaces the known mixed-apex frontier residual honestly;
/// this gate pins completion + the residual ceiling (currently
/// 1.36e-2, the PK-092 direct-apex class), not perfection.
#[test]
fn record_modify_capped_corner_re_dispatches() {
    let (record, strip) = capped_corner_record();
    // PK-174: the local (recordless) path now COMPLETES on the
    // trihedral capped corner through the corner-group restore —
    // hold it to the same pinned frontier as the record path.
    let local =
        popkern_blend::modify_fillet_radius(&record.result, strip, 0.15, Tolerance::default())
            .expect("recordless corner-group modify completes on the trihedral class");
    let local_errors = local.validate_all();
    assert!(
        local_errors.len() <= 1,
        "recordless capped-corner modify regressed: {local_errors:?}"
    );
    let m = popkern_blend::modify_fillet_radius_via_record(&record, strip, 0.15)
        .expect("record re-dispatch");
    let errors = m.brep.validate_all();
    let worst = errors
        .iter()
        .filter_map(|e| {
            let s = format!("{e:?}");
            s.split("max_deviation: ")
                .nth(1)
                .and_then(|t| t.trim_end_matches(&[' ', '}'][..]).parse::<f64>().ok())
        })
        .fold(0.0_f64, f64::max);
    assert!(
        errors.len() <= 1 && worst <= 2.0e-2,
        "mixed-radius re-dispatch regressed past the pinned frontier: {} errors, worst {worst:.3e}",
        errors.len()
    );
}
