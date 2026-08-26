//! PK-326 row 477: face collapse. Features thinner than the offset
//! dissolve out of the skin through Euler operators (`kef` merges,
//! `kev` spur retirement, terminal `kvr`) with a `FaceCollapsed`
//! warning naming the vanished faces and the measured mechanism; the
//! surviving skin is the true Minkowski erosion. Folded geometry never
//! serves: anything the dissolution cannot resolve still declines
//! through the chord-inversion gate.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::thin_rib;
use popkern_mesh::offset_oracle::{
    OffsetOracleOptions, OffsetSemantics, classify_offset_agreement,
};
use popkern_mesh::point_query::PointMeshIndex;
use popkern_mesh::{TessellationOptions, tessellate};
use popkern_shell::offset::offset_solid;
use popkern_shell::thicken::thicken;
use popkern_shell::{ShellError, ShellOptions, ShellWarning};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn index_of(brep: &BRep) -> PointMeshIndex {
    let mesh = tessellate(brep, &TessellationOptions::default()).expect("tessellate");
    PointMeshIndex::build(&mesh)
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

fn collapsed_faces(warnings: &[ShellWarning]) -> Option<&Vec<popkern_core::topology::FaceHandle>> {
    warnings.iter().find_map(|w| match w {
        ShellWarning::FaceCollapsed { faces, .. } => Some(faces),
        _ => None,
    })
}

/// Erosion past half the rib width (0.03) dissolves the whole rib:
/// the result is EXACTLY the eroded plate box, with the five rib
/// faces named in the warning and the erosion predicate holding with
/// zero exclusions.
#[test]
fn thin_rib_erosion_collapses_the_rib_exactly() {
    let tr = thin_rib(tol()).unwrap().brep;
    let d = 0.04;
    let result = offset_solid(&tr, -d, &ShellOptions::default(), tol()).expect("erosion");
    assert!(result.brep.validate_all().is_empty());

    let faces = collapsed_faces(&result.warnings).expect("FaceCollapsed warning");
    assert_eq!(faces.len(), 5, "the five rib faces vanish: {faces:?}");

    // The eroded body is a plain box.
    assert_eq!(result.brep.faces.len(), 6);
    assert_eq!(result.brep.edges.len(), 12);
    assert_eq!(result.brep.vertices.len(), 8);

    // Exact closed form: [d,2−d] × [d,1−d] × [d,0.2−d].
    let expected = (2.0 - 2.0 * d) * (1.0 - 2.0 * d) * (0.2 - 2.0 * d);
    let measured = mesh_volume(&result.brep);
    assert!(
        (measured - expected).abs() / expected < 1e-9,
        "volume {measured} vs {expected}"
    );

    let report = classify_offset_agreement(
        &index_of(&tr),
        &index_of(&result.brep),
        OffsetSemantics::Shrink { distance: d },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert_eq!(report.excluded, 0);
    assert!(
        report.is_clean(),
        "erosion oracle disagreements: {} (first {:?})",
        report.disagreement_count,
        report.disagreements.first()
    );
}

/// Thicken past the rib width: the wall solid keeps the rib entirely
/// (it is thinner than the wall) while the inner void is the eroded
/// box. Serves with the warning, validates, and matches the Wall
/// predicate with zero exclusions.
#[test]
fn thin_rib_thicken_past_rib_width_serves_wall() {
    let tr = thin_rib(tol()).unwrap().brep;
    let t = 0.08; // per side t/2 = 0.04 > rib half-width 0.03
    let result = thicken(&tr, t, &ShellOptions::default(), tol()).expect("thicken");
    assert!(result.brep.validate_all().is_empty());
    assert!(
        collapsed_faces(&result.warnings).is_some(),
        "warnings: {:?}",
        result.warnings
    );

    let report = classify_offset_agreement(
        &index_of(&tr),
        &index_of(&result.brep),
        OffsetSemantics::Wall {
            inward: t / 2.0,
            outward: t / 2.0,
        },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert_eq!(report.excluded, 0);
    assert!(
        report.is_clean(),
        "wall oracle disagreements: {} (first {:?})",
        report.disagreement_count,
        report.disagreements.first()
    );
}

/// Below the threshold the rib survives (no warning, 15 faces: the
/// footprint ring rounds via the mixed-corner arc joins); exactly AT
/// the degenerate threshold the fold declines — never serves.
#[test]
fn collapse_threshold_is_honest() {
    let tr = thin_rib(tol()).unwrap().brep;

    let under = offset_solid(&tr, -0.02, &ShellOptions::default(), tol()).expect("serves");
    assert!(collapsed_faces(&under.warnings).is_none());
    assert_eq!(under.brep.faces.len(), 15, "rib survives with 4 gap strips");
    assert!(under.brep.validate_all().is_empty());

    let at = offset_solid(&tr, -0.03, &ShellOptions::default(), tol());
    assert!(
        matches!(at, Err(ShellError::OffsetTooLarge { .. })),
        "degenerate zero-width rib must decline, got {at:?}"
    );
}

/// Whole-body over-erosion has no surviving skin to dissolve into:
/// the chord-inversion gate still declines with OffsetTooLarge.
#[test]
fn whole_body_collapse_still_declines() {
    let cube = BRep::unit_cube(tol()).unwrap();
    for d in [0.5, 0.6] {
        let r = offset_solid(&cube, -d, &ShellOptions::default(), tol());
        assert!(
            matches!(r, Err(ShellError::OffsetTooLarge { .. })),
            "cube -{d}: {r:?}"
        );
    }
}
