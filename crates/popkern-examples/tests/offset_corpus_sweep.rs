//! Corpus sweep for the offset campaign (PK-323, design §9): every
//! fixture builds, carries its derived topology/census pins, validates
//! fully, and measures its analytic volume. The campaign-level gate
//! lives here too: each closed fixture either SERVES an inward thicken
//! (validator + SDF oracle clean) or DECLINES with a structured
//! diagnostic — zero silent wrong geometry. Serve/decline expectations
//! are pinned per fixture; a later campaign stage that widens the
//! serving envelope flips its pin deliberately, with oracle evidence.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::{corpus, surface_census};
use popkern_mesh::offset_oracle::{
    OffsetOracleOptions, OffsetSemantics, classify_offset_agreement,
};
use popkern_mesh::point_query::PointMeshIndex;
use popkern_mesh::{TessellationOptions, tessellate};
use popkern_shell::{ShellOptions, thicken::thicken};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn mesh_volume(brep: &BRep) -> f64 {
    let mesh = tessellate(brep, &TessellationOptions::default()).expect("tessellate");
    let mut vol = 0.0_f64;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (vol / 6.0).abs()
}

#[test]
fn corpus_topology_census_validation_and_volume_pins() {
    for f in corpus(tol()).expect("corpus builds") {
        assert_eq!(f.brep.faces.len(), f.expected_faces, "{}: faces", f.name);
        assert_eq!(f.brep.edges.len(), f.expected_edges, "{}: edges", f.name);
        assert_eq!(
            f.brep.vertices.len(),
            f.expected_vertices,
            "{}: vertices",
            f.name
        );
        assert_eq!(
            surface_census(&f.brep),
            f.expected_census,
            "{}: surface census",
            f.name
        );
        let errors = f.brep.validate_all();
        assert!(errors.is_empty(), "{}: validate_all: {errors:?}", f.name);
        if let Some(expected) = f.analytic_volume {
            let measured = mesh_volume(&f.brep);
            let rel = (measured - expected).abs() / expected;
            assert!(
                rel <= f.volume_band,
                "{}: volume {measured:.6} vs analytic {expected:.6} (rel {rel:.5}, band {})",
                f.name,
                f.volume_band
            );
        }
    }
}

/// Which fixtures the CURRENT serving envelope thickens. Every other
/// fixture must decline with a structured error. Widening stages flip
/// entries to `true` alongside their oracle evidence (PK-324 rows).
fn expected_to_serve(name: &str) -> bool {
    // PK-324 flips: filleted_cube (parallel-map arm serves its 48
    // G1-tangential edges), bracket_bore_countersink (coaxial-profile
    // arm serves the bore×countersink cylinder×cone junction circle),
    // torus_lens (parallel-map arm serves its same-torus patch-border
    // circles). sphere_capped_rod: the parallel-map arm gained
    // endpoint re-trim (its sphere-meridian seam pieces end on the
    // transversal junction circle, whose vertex slides along the
    // meridian under offset).
    //
    // PK-325 flip: wavy_plate — the Stage-3 acceptance body. Its
    // convex top-boundary edges separate under the outward offset;
    // the arc-join stage now closes them with rolling-ball cylinder
    // strips (the top boundary rows are straight) and sphere octant
    // corner caps, verified below by validate_all + the SDF oracle
    // on the served body.
    matches!(
        name,
        "l_bracket"
            | "thin_rib"
            | "high_valence_pyramid"
            | "filleted_cube"
            | "bracket_bore_countersink"
            | "torus_lens"
            | "sphere_capped_rod"
            | "wavy_plate"
    )
}

#[test]
fn corpus_campaign_gate_thicken_serves_clean_or_declines_structured() {
    let t = 0.02; // below the thin rib's collapse threshold (0.03)
    for f in corpus(tol()).expect("corpus builds") {
        if !f.watertight {
            // Open sheets thicken through the free-rim wall path
            // (PK-324): serve a CLOSED solid, validate_all-clean, with
            // the exact slab volume. The SDF Wall oracle needs a
            // watertight SOURCE for its containment side, so the open
            // sheet is gated on volume + closure instead.
            let result = thicken(&f.brep, t, &ShellOptions::default(), tol())
                .unwrap_or_else(|e| panic!("{}: open-sheet thicken declined: {e}", f.name));
            let errors = result.brep.validate_all();
            assert!(errors.is_empty(), "{}: served invalid: {errors:?}", f.name);
            let free = result
                .brep
                .half_edges
                .iter()
                .filter(|(_, he)| !result.brep.loops.contains_key(he.loop_))
                .count();
            assert_eq!(free, 0, "{}: free rims remain", f.name);
            // Closed-form (PK-326 rim corners): the OUTER skin serves
            // rounded — full Minkowski dilation of the box by t/2
            // intersected with the opening half-space z ≤ 1:
            //   V_outer = 1 + 5t' + 2π t'² + (2π/3) t'³ , t' = t/2
            // (8 quarter-round edges, 4 sphere-octant corners, rim
            // arcs at the opening plane). The inner void stays the
            // sharp box (1−t)²(1−t/2): erosion keeps convex edges
            // sharp. Rim walls remain flat in the opening plane.
            let tp = t / 2.0;
            let outer_vol = 1.0
                + 5.0 * tp
                + 2.0 * std::f64::consts::PI * tp * tp
                + (2.0 * std::f64::consts::PI / 3.0) * tp * tp * tp;
            let expected = outer_vol - (1.0 - t) * (1.0 - t) * (1.0 - t / 2.0);
            let measured = mesh_volume(&result.brep);
            let rel = (measured - expected).abs() / expected;
            assert!(
                rel < 5e-3,
                "{}: wall volume {measured:.6} vs analytic {expected:.6} (rel {rel:.4})",
                f.name
            );
            assert!(
                !result
                    .warnings
                    .iter()
                    .any(|w| matches!(w, popkern_shell::ShellWarning::JoinStyleDowngraded { .. })),
                "{}: rim corners must serve without downgrade: {:?}",
                f.name,
                result.warnings
            );
            println!("[gate] {} SERVES (open sheet): volume rel {rel:.5}", f.name);
            continue;
        }
        match thicken(&f.brep, t, &ShellOptions::default(), tol()) {
            Ok(result) => {
                assert!(
                    expected_to_serve(f.name),
                    "{}: thicken now serves but was pinned declining — verify with the \
                     oracle and flip the pin deliberately",
                    f.name
                );
                let errors = result.brep.validate_all();
                assert!(errors.is_empty(), "{}: served invalid: {errors:?}", f.name);
                let source_mesh =
                    tessellate(&f.brep, &TessellationOptions::default()).expect("tessellate");
                let result_mesh =
                    tessellate(&result.brep, &TessellationOptions::default()).expect("tessellate");
                let source = PointMeshIndex::build(&source_mesh);
                let served = PointMeshIndex::build(&result_mesh);
                // PK-326: no exclusion zones. The old reflex-wedge
                // exclusions covered the SHARP inner corners the
                // downgraded (Intersection) joins kept near concave
                // edges; mixed corners now serve under Arc, so the
                // rounded joins are distance-faithful everywhere.
                let report = classify_offset_agreement(
                    &source,
                    &served,
                    OffsetSemantics::Wall {
                        inward: t / 2.0,
                        outward: t / 2.0,
                    },
                    &OffsetOracleOptions::default(),
                    |_| false,
                );
                // PK-326 warning-retirement pin: the E1 fixtures serve
                // with no style downgrades.
                if matches!(f.name, "l_bracket" | "thin_rib") {
                    assert!(
                        !result.warnings.iter().any(|w| matches!(
                            w,
                            popkern_shell::ShellWarning::JoinStyleDowngraded { .. }
                                | popkern_shell::ShellWarning::JoinStyleDowngradedAtEdge { .. }
                        )),
                        "{}: unexpected join downgrade: {:?}",
                        f.name,
                        result.warnings
                    );
                }
                assert!(
                    report.is_clean(),
                    "{}: SDF oracle disagreements: {} (first: {:?})",
                    f.name,
                    report.disagreement_count,
                    report.disagreements.first()
                );
                assert!(
                    report.classified_fraction() > 0.5,
                    "{}: oracle classified too little: {report:?}",
                    f.name
                );
                assert!(
                    (report.excluded as f64) <= 0.05 * report.total_samples as f64,
                    "{}: exclusion zone too large: {} of {}",
                    f.name,
                    report.excluded,
                    report.total_samples
                );
                println!(
                    "[gate] {} SERVES: oracle {} agreements, {} boundary-skipped, {} excluded",
                    f.name, report.agreements, report.boundary_skipped, report.excluded
                );
            }
            Err(e) => {
                assert!(
                    !expected_to_serve(f.name),
                    "{}: pinned serving but declined: {e}",
                    f.name
                );
                println!("[gate] {} declines (structured): {e}", f.name);
            }
        }
    }
}
