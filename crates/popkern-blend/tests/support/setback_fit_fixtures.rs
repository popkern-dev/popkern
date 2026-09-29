//! Shared fixtures and assertions for the variational corner fit tests
//! (PK-574): fillet helpers, the fit-record accessors and the row-919
//! serve assertions.
#![allow(dead_code)]

use popkern_blend::dispatch_overrides::with_kp20_cap_forced_to_decline;
use popkern_blend::{
    ApexConstructionDiagnostic, BlendResult, ContinuityLevel, CornerMethod, FilletOptions,
    FilletRadius, FilletSpec, SetbackFitDiagnostic, SetbackFitMode, SetbackFitOutcome,
    fillet_edges,
};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

pub fn edges_at(body: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
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

pub fn fillet(
    body: &BRep,
    vertex: VertexHandle,
    radii: &[f64],
    continuity: ContinuityLevel,
) -> BlendResult {
    try_fillet(body, vertex, radii, continuity).expect("fillet")
}

pub fn try_fillet(
    body: &BRep,
    vertex: VertexHandle,
    radii: &[f64],
    continuity: ContinuityLevel,
) -> Result<BlendResult, String> {
    let specs = edges_at(body, vertex)
        .into_iter()
        .zip(radii.iter().cycle())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect::<Vec<_>>();
    let options = FilletOptions {
        continuity,
        ..Default::default()
    };
    fillet_edges(body, &specs, &options, Tolerance::default(), None).map_err(|e| {
        let text = format!("{e:?}");
        text[..text.len().min(300)].to_string()
    })
}

pub fn cube_corner(radii: &[f64]) -> BlendResult {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(1.0, 1.0, 1.0)) < 1e-9)
        .map(|(h, _)| h)
        .expect("corner");
    fillet(&cube, corner, radii, ContinuityLevel::G2)
}

pub fn fit_records(result: &BlendResult) -> Vec<&SetbackFitDiagnostic> {
    result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            ApexConstructionDiagnostic::SetbackFit(r) => Some(r),
            _ => None,
        })
        .collect()
}

/// Commercial gates (`ContinuityGates::commercial`).
pub const G0_GATE: f64 = 1e-5;
pub const G1_GATE: f64 = 5e-2;
pub const G2_GATE: f64 = 5e-2;

/// The single fit record of a one-corner fillet.
pub fn only_record(name: &str, result: &BlendResult) -> SetbackFitDiagnostic {
    let records = fit_records(result);
    assert_eq!(records.len(), 1, "{name}: one fit record, got {records:#?}");
    records[0].clone()
}

/// Row 919: the fit serves at `level` in `mode`, its certificate
/// covers every cap seam at that level, the spokes are G2 (exact by
/// construction), and the committed body validates with the corner's
/// sectors recorded as fitted plate patches.
pub fn assert_fit_serves(
    name: &str,
    result: &BlendResult,
    mode: SetbackFitMode,
    level: ContinuityLevel,
) {
    let record = only_record(name, result);
    assert_record_serves(name, &record, mode, level);
    let fitted = result
        .corner_patches
        .iter()
        .filter(|p| p.method == CornerMethod::PlateMin)
        .inspect(|p| assert_eq!(p.achieved_continuity, level, "{name}: {p:?}"))
        .count();
    assert!(fitted >= 3, "{name}: {fitted} fitted sectors");
    // Row 880: stable ownership — one fitted sector patch per sector
    // index, each an ordinary NURBS face.
    let mut sectors: Vec<usize> = result
        .corner_patches
        .iter()
        .filter(|p| p.method == CornerMethod::PlateMin)
        .map(|p| {
            let prov = p.sub_patch_provenance.as_ref().expect("sector provenance");
            assert_eq!(
                prov.family,
                popkern_blend::types::CornerSubPatchFamily::VariationalFit,
                "{name}"
            );
            assert_eq!(
                prov.role,
                popkern_blend::types::CornerSubPatchRole::Sector,
                "{name}"
            );
            assert_eq!(prov.sector_count, fitted, "{name}");
            assert!(
                matches!(
                    result.brep.faces[p.face]
                        .surface
                        .as_ref()
                        .map(|s| &s.geometry),
                    Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_))
                ),
                "{name}: fitted sector is a NURBS face"
            );
            prov.sector_index
        })
        .collect();
    sectors.sort_unstable();
    assert_eq!(
        sectors,
        (0..fitted).collect::<Vec<_>>(),
        "{name}: sector indices"
    );
    result
        .brep
        .validate()
        .unwrap_or_else(|e| panic!("{name}: {e}"));
}

/// One corner's record: served at `level` by `mode`, with the measured
/// bounds row 919 pins.
pub fn assert_record_serves(
    name: &str,
    record: &SetbackFitDiagnostic,
    mode: SetbackFitMode,
    level: ContinuityLevel,
) {
    assert_eq!(record.served, Some(level), "{name}: {record:#?}");
    let served = record.attempts.last().expect("an attempt");
    assert_eq!(served.mode, mode, "{name}: {record:#?}");
    let SetbackFitOutcome::Served(c) = &served.outcome else {
        panic!("{name}: last attempt not served: {record:#?}");
    };
    assert_eq!(c.level, level, "{name}");
    // E6: served caps measured at regularity 0.131 and above; the core
    // fold gate declines below 0.02.
    assert!(
        c.min_regularity >= 0.1,
        "{name}: regularity {}",
        c.min_regularity
    );
    // Measured bounds (refined past the gates): rims G1 more than three
    // decades inside the gate. Spokes are G2 by construction (every
    // certified lower witness is 0); their upper bounds are enclosure
    // width at the interval budget (worst measured 1.5e-2, the wide
    // mixed cube), so they are held to the gate itself.
    assert!(
        c.rim[0] <= G0_GATE && c.rim[1] <= 1e-5,
        "{name}: rim {:?}",
        c.rim
    );
    assert!(
        c.spoke[0] <= G0_GATE && c.spoke[1] <= 1e-5 && c.spoke[2] <= G2_GATE,
        "{name}: spoke {:?}",
        c.spoke
    );
    if level == ContinuityLevel::G2 {
        assert!(c.rim[2] <= G2_GATE, "{name}: rim {:?}", c.rim);
    }
}

pub fn pyramid_apex(
    n: usize,
    height: f64,
    radii: &[f64],
    continuity: ContinuityLevel,
) -> BlendResult {
    let body = BRep::n_gonal_pyramid(n, 1.0, height, Tolerance::default()).expect("pyramid");
    let apex = body
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex");
    fillet(&body, apex, radii, continuity)
}

/// A named fixture with the mode and level the fit must serve at.
pub type Fixture = (
    &'static str,
    Box<dyn Fn() -> BlendResult>,
    SetbackFitMode,
    ContinuityLevel,
);

/// Row 876: every setback corner the KP20 cap serves is served by the
/// fit when KP20 is forced to decline at its cap stage — at certified
/// G1 by the tangent fit, or at certified G2 by the curvature-weighted
/// fit where the corner admits it (the setback corner theorem, design
/// doc: G2 strips with nonzero curvature slope admit no C^3 G2 cap).
pub fn census(fixtures: Vec<Fixture>) {
    for (name, build, mode, level) in fixtures {
        let result = with_kp20_cap_forced_to_decline(build);
        assert_fit_serves(name, &result, mode, level);
    }
}

/// A unit cube with its (+,+,+) corner notched out, and the saddle
/// vertex at (0.6, 0.6, 0.6) (the PK-527 fixture).
pub fn notched_saddle() -> (BRep, VertexHandle) {
    use popkern_core::geometry::transform::Transform3;
    let tol = Tolerance::default();
    let notch = Transform3::translation(0.6, 0.6, 0.6).compose(&Transform3::scale_uniform(0.6));
    let tool = BRep::unit_cube(tol)
        .expect("cube")
        .transform(&notch)
        .expect("notch");
    let body = popkern_bool::subtraction(BRep::unit_cube(tol).expect("cube"), tool, tol, None)
        .expect("notched cube");
    let saddle = body
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(0.6, 0.6, 0.6)) < 1e-6)
        .map(|(h, _)| h)
        .expect("saddle vertex");
    (body, saddle)
}
