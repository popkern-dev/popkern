//! PK-326 row 476 (exploration E1): mixed convex/concave corners serve
//! under `JoinStyle::Arc` via the two-sector transversal peel — the
//! measured probe refuted both pre-forked mechanisms (no sphere patch
//! exists: the polar cone has empty interior), and the connector edge
//! is the transversal intersection of the two cut-wedge surfaces.
//!
//! The gates run at d = 0.1, where a sharp (extended-intersection)
//! join deviates from the Minkowski predicate by d·(√2−1) ≈ 0.041 —
//! twice the oracle's 0.02 skip band — so agreement here genuinely
//! pins the rounded transversal geometry, not band arithmetic.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::{l_bracket, surface_census, thin_rib};
use popkern_mesh::offset_oracle::{
    OffsetOracleOptions, OffsetSemantics, classify_offset_agreement,
};
use popkern_mesh::point_query::PointMeshIndex;
use popkern_mesh::{TessellationOptions, tessellate};
use popkern_shell::ShellOptions;
use popkern_shell::offset::offset_solid;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn index_of(brep: &BRep) -> PointMeshIndex {
    let mesh = tessellate(brep, &TessellationOptions::default()).expect("tessellate");
    PointMeshIndex::build(&mesh)
}

fn assert_clean(report: &popkern_mesh::offset_oracle::OffsetOracleReport, label: &str) {
    assert!(
        report.is_clean(),
        "{label}: SDF oracle disagreements: {} (first: {:?})",
        report.disagreement_count,
        report.disagreements.first()
    );
    assert!(
        report.classified_fraction() > 0.5,
        "{label}: oracle classified too little: {report:?}"
    );
}

/// Saddle corners (2 gap + 1 resolvable edge): the L-bracket outward
/// offset serves Minkowski-exact with zero warnings and zero oracle
/// exclusions. Topology: 17 strips, 10 sphere-triangle corners, 2
/// saddle corners with ellipse connectors (strip×strip transversal).
#[test]
fn l_bracket_outward_offset_serves_minkowski_exact() {
    let lb = l_bracket(tol()).unwrap().brep;
    let d = 0.1;
    let result = offset_solid(&lb, d, &ShellOptions::default(), tol()).expect("arc offset");
    assert!(
        result.warnings.is_empty(),
        "expected zero containment warnings, got {:?}",
        result.warnings
    );
    assert!(result.brep.validate_all().is_empty());

    // Topology pins: 8 faces + 17 strips + 10 caps; per gap edge +1E,
    // per trihedral corner +2V+3E, per saddle corner +1V+1E.
    assert_eq!(result.brep.faces.len(), 35, "faces");
    assert_eq!(result.brep.edges.len(), 67, "edges");
    assert_eq!(result.brep.vertices.len(), 34, "vertices");
    assert_eq!(
        surface_census(&result.brep),
        vec![("cylinder", 17), ("plane", 8), ("sphere", 10)],
        "surface census"
    );

    let report = classify_offset_agreement(
        &index_of(&lb),
        &index_of(&result.brep),
        OffsetSemantics::Grow { distance: d },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert_eq!(report.excluded, 0);
    assert_clean(&report, "l_bracket outward (Grow 0.1)");
}

/// T-corners (1 gap + 2 resolvable edges): the L-bracket inward offset
/// erodes with rounded concave joins — the reflex edge's strip
/// dead-ends into the cap offset planes along exact circle connectors.
/// A sharp join misses the rounded wedge material by ≈ 0.041 here, so
/// zero-exclusion Shrink agreement is unreachable without the arc.
#[test]
fn l_bracket_inward_offset_serves_erosion_exact() {
    let lb = l_bracket(tol()).unwrap().brep;
    let d = 0.1;
    let result = offset_solid(&lb, -d, &ShellOptions::default(), tol()).expect("arc erosion");
    assert!(
        result.warnings.is_empty(),
        "expected zero containment warnings, got {:?}",
        result.warnings
    );
    assert!(result.brep.validate_all().is_empty());

    // One gap edge (the reflex vertical), two T-corners: +1E for the
    // strip, +1V+1E per T-corner.
    assert_eq!(result.brep.faces.len(), 9, "faces");
    assert_eq!(result.brep.edges.len(), 21, "edges");
    assert_eq!(result.brep.vertices.len(), 14, "vertices");
    assert_eq!(
        surface_census(&result.brep),
        vec![("cylinder", 1), ("plane", 8)],
        "surface census"
    );

    let report = classify_offset_agreement(
        &index_of(&lb),
        &index_of(&result.brep),
        OffsetSemantics::Shrink { distance: d },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert_eq!(report.excluded, 0);
    assert_clean(&report, "l_bracket inward (Shrink 0.1)");
}

/// Dangling-gap T-corners at the rib base: the rib's vertical strip
/// corners dead-end into the plate-top offset plane along quarter
/// circles; the footprint stays resolvable. Zero warnings, exact Grow.
#[test]
fn thin_rib_outward_offset_serves_minkowski_exact() {
    let tr = thin_rib(tol()).unwrap().brep;
    let d = 0.1;
    let result = offset_solid(&tr, d, &ShellOptions::default(), tol()).expect("arc offset");
    assert!(
        result.warnings.is_empty(),
        "expected zero containment warnings, got {:?}",
        result.warnings
    );
    assert!(result.brep.validate_all().is_empty());

    // 20 gap edges (plate ring 12, rib verticals 4, rib top 4), 12
    // trihedral corners, 4 T-corners at the rib base.
    assert_eq!(result.brep.faces.len(), 43, "faces");
    assert_eq!(result.brep.edges.len(), 84, "edges");
    assert_eq!(result.brep.vertices.len(), 44, "vertices");
    assert_eq!(
        surface_census(&result.brep),
        vec![("cylinder", 20), ("plane", 11), ("sphere", 12)],
        "surface census"
    );

    let report = classify_offset_agreement(
        &index_of(&tr),
        &index_of(&result.brep),
        OffsetSemantics::Grow { distance: d },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert_eq!(report.excluded, 0);
    assert_clean(&report, "thin_rib outward (Grow 0.1)");
}

/// Determinism: the mixed-corner construction is bit-identical across
/// repeated runs (sorted surgery order, no ambient state).
#[test]
fn mixed_corner_construction_is_bit_identical() {
    let lb = l_bracket(tol()).unwrap().brep;
    let run = || {
        let result = offset_solid(&lb, 0.1, &ShellOptions::default(), tol()).expect("offset");
        let mut bits: Vec<[u64; 3]> = result
            .brep
            .vertices
            .iter()
            .map(|(_, v)| {
                [
                    v.position.x.to_bits(),
                    v.position.y.to_bits(),
                    v.position.z.to_bits(),
                ]
            })
            .collect();
        bits.sort();
        bits
    };
    let first = run();
    for _ in 0..2 {
        assert_eq!(run(), first, "vertex bits differ between runs");
    }
}

/// Dilated cone apex (PK-326): the Minkowski dilation rounds the apex
/// with the sphere of radius d centered at the ORIGINAL apex, tangent
/// to the offset cone along an exact circle (radius d·cosβ). The
/// `ConeApexSharp` warning retires; erosion keeps the sharp apex
/// (exact) with no cap and no warning.
#[test]
fn cone_apex_serves_sphere_cap_under_dilation() {
    let cone = BRep::cone(1.0, 1.5, tol()).unwrap();
    let d = 0.1;
    let result = offset_solid(&cone, d, &ShellOptions::default(), tol()).expect("arc offset");
    assert!(
        result.warnings.is_empty(),
        "apex warning must retire: {:?}",
        result.warnings
    );
    assert!(result.brep.validate_all().is_empty());
    assert_eq!(
        surface_census(&result.brep),
        vec![("cone", 1), ("plane", 1), ("sphere", 1), ("torus", 2)],
        "lateral + base + apex cap + split rim torus window"
    );

    let report = classify_offset_agreement(
        &index_of(&cone),
        &index_of(&result.brep),
        OffsetSemantics::Grow { distance: d },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert_eq!(report.excluded, 0);
    assert_clean(&report, "cone outward (Grow 0.1, apex sphere cap)");

    // Erosion: apex stays sharp and exact, no cap, no warning.
    let eroded = offset_solid(&cone, -d, &ShellOptions::default(), tol()).expect("erosion");
    assert!(eroded.warnings.is_empty(), "{:?}", eroded.warnings);
    assert!(eroded.brep.validate_all().is_empty());
    assert_eq!(
        surface_census(&eroded.brep),
        vec![("cone", 1), ("plane", 1)],
        "eroded cone keeps its sharp apex"
    );
    let report = classify_offset_agreement(
        &index_of(&cone),
        &index_of(&eroded.brep),
        OffsetSemantics::Shrink { distance: d },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert_eq!(report.excluded, 0);
    assert_clean(&report, "cone inward (Shrink 0.1)");
}

/// Rim (free-boundary) corners (PK-326): thickening an open sheet
/// serves the OUTER skin rounded — strips on every convex edge,
/// sphere triangles at closed corners, and canal end-circle arcs at
/// the four boundary corners — with flat rim walls (quads plus corner
/// sectors) closing the solid. No downgrade warnings.
#[test]
fn open_sheet_thicken_serves_rounded_with_rim_arc_corners() {
    use popkern_examples::offset_corpus::open_sheet;
    use popkern_shell::thicken::thicken;

    let sheet = open_sheet(tol()).unwrap().brep;
    let t = 0.05;
    let result = thicken(&sheet, t, &ShellOptions::default(), tol()).expect("thicken");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| matches!(w, popkern_shell::ShellWarning::JoinStyleDowngraded { .. })),
        "rim corners must serve without downgrade: {:?}",
        result.warnings
    );
    assert!(result.brep.validate_all().is_empty());
    // 17 outer (5 faces + 8 strips + 4 sphere triangles) + 5 inner +
    // 4 rim wall quads + 4 rim corner sectors.
    assert_eq!(result.brep.faces.len(), 30, "faces");
    // Closed solid: no loopless half-edges survive.
    let free = result
        .brep
        .half_edges
        .iter()
        .filter(|(_, he)| !result.brep.loops.contains_key(he.loop_))
        .count();
    assert_eq!(free, 0, "free half-edges remain");
    // Exact closed form: full Minkowski dilation of the box by t/2
    // intersected with the opening half-space, minus the sharp inner
    // void.
    let tp = t / 2.0;
    let outer_vol = 1.0
        + 5.0 * tp
        + 2.0 * std::f64::consts::PI * tp * tp
        + (2.0 * std::f64::consts::PI / 3.0) * tp * tp * tp;
    let expected = outer_vol - (1.0 - t) * (1.0 - t) * (1.0 - t / 2.0);
    let mesh = tessellate(&result.brep, &TessellationOptions::default()).expect("tessellate");
    let mut vol = 0.0_f64;
    for tri in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[tri[0]].position,
            mesh.vertices[tri[1]].position,
            mesh.vertices[tri[2]].position,
        );
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    let measured = (vol / 6.0).abs();
    let rel = (measured - expected).abs() / expected;
    assert!(
        rel < 5e-3,
        "wall volume {measured:.6} vs rounded closed form {expected:.6} (rel {rel:.4})"
    );
}
