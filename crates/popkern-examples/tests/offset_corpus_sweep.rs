//! Corpus sweep for the offset campaign (PK-323, design §9): every
//! fixture builds, carries its derived topology/census pins, validates
//! fully, and measures its analytic volume. The campaign-level gate
//! lives here too: each closed fixture either SERVES an inward thicken
//! (validator + SDF oracle clean) or DECLINES with a structured
//! diagnostic — zero silent wrong geometry. Serve/decline expectations
//! are pinned per fixture; a later campaign stage that widens the
//! serving envelope flips its pin deliberately, with oracle evidence.

use popkern_core::geometry::Point3;
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

/// Reflex-wedge exclusion for the sharp (intersection-join) inward
/// thicken: near a concave edge the sharp inner corner legitimately
/// keeps material the pure distance band would not (reach `t·√2` for
/// right angles). Points within `radius` of any listed segment are
/// excluded — counted, never silent.
fn near_segments(p: Point3, segments: &[(Point3, Point3)], radius: f64) -> bool {
    segments.iter().any(|(a, b)| {
        let ab = *b - *a;
        let len_sq = ab.dot(&ab);
        let t = if len_sq > 0.0 {
            ((p - *a).dot(&ab) / len_sq).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let q = *a + ab * t;
        let d = p - q;
        d.dot(&d).sqrt() < radius
    })
}

/// Concave-edge segments per fixture (for the exclusion closure).
fn concave_edges(name: &str) -> Vec<(Point3, Point3)> {
    match name {
        "l_bracket" => vec![(Point3::new(0.4, 0.4, 0.0), Point3::new(0.4, 0.4, 0.6))],
        "thin_rib" => {
            // Rib footprint ring on the plate top (z = 0.2).
            let (x0, x1, y0, y1, z) = (0.2, 1.8, 0.47, 0.53, 0.2);
            vec![
                (Point3::new(x0, y0, z), Point3::new(x1, y0, z)),
                (Point3::new(x1, y0, z), Point3::new(x1, y1, z)),
                (Point3::new(x1, y1, z), Point3::new(x0, y1, z)),
                (Point3::new(x0, y1, z), Point3::new(x0, y0, z)),
            ]
        }
        _ => vec![],
    }
}

/// Which fixtures the CURRENT serving envelope thickens (sharp inward
/// joins over exact analytic offsets). Every other fixture must
/// decline with a structured error. Widening stages flip entries to
/// `true` alongside their oracle evidence (PK-324 rows).
fn expected_to_serve(name: &str) -> bool {
    // PK-324 flips: filleted_cube (parallel-map arm serves its 48
    // G1-tangential edges), bracket_bore_countersink (coaxial-profile
    // arm serves the bore×countersink cylinder×cone junction circle),
    // torus_lens (parallel-map arm serves its same-torus patch-border
    // circles). wavy_plate still declines — measured: the INWARD
    // ladder recomputes fully (general SSI arm), but the OUTWARD
    // offset of its convex top-boundary edges separates the surfaces
    // and the fitted NURBS offset is domain-bounded, so no
    // topology-preserving intersection exists; closing that gap is
    // the Stage-3 join/extension machinery (PK-325's acceptance body).
    // sphere_capped_rod: the parallel-map arm gained endpoint re-trim
    // (its sphere-meridian seam pieces end on the transversal junction
    // circle, whose vertex slides along the meridian under offset).
    matches!(
        name,
        "l_bracket"
            | "thin_rib"
            | "high_valence_pyramid"
            | "filleted_cube"
            | "bracket_bore_countersink"
            | "torus_lens"
            | "sphere_capped_rod"
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
            // Closed-form: outer box (1+t)²(1+t/2) minus inner void
            // (1−t)²(1−t/2) — the open [0,1]³ box's five slabs with
            // sharp rim walls at the opening plane z = 1.
            let expected =
                (1.0 + t) * (1.0 + t) * (1.0 + t / 2.0) - (1.0 - t) * (1.0 - t) * (1.0 - t / 2.0);
            let measured = mesh_volume(&result.brep);
            let rel = (measured - expected).abs() / expected;
            assert!(
                rel < 0.02,
                "{}: wall volume {measured:.6} vs analytic {expected:.6} (rel {rel:.4})",
                f.name
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
                let wedges = concave_edges(f.name);
                let report = classify_offset_agreement(
                    &source,
                    &served,
                    OffsetSemantics::Wall {
                        inward: t / 2.0,
                        outward: t / 2.0,
                    },
                    &OffsetOracleOptions::default(),
                    |p| near_segments(p, &wedges, 2.5 * t),
                );
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
