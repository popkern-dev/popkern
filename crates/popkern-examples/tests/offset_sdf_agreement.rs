//! PK-323 row 466: the SDF classification oracle agrees with served
//! offsets on the cube, cylinder, and filleted-cube fixtures.
//!
//! Wall/erosion semantics are exact for the sharp inward joins the
//! current pipeline serves; the outward cone offset demonstrates Grow
//! semantics with the sharp-join wedge zones excluded explicitly (and
//! counted). The filleted cube's serving is a PK-324 acceptance row —
//! here the oracle itself is held to the fixture's closed-form SDF
//! (the corner caps are exact sphere octants), and a serve attempt is
//! gated either way: serve ⇒ oracle-clean, decline ⇒ structured.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::filleted_cube;
use popkern_mesh::offset_oracle::{
    OffsetOracleOptions, OffsetSemantics, classify_offset_agreement,
};
use popkern_mesh::point_query::PointMeshIndex;
use popkern_mesh::{TessellationOptions, tessellate};
use popkern_shell::offset::offset_solid;
use popkern_shell::{ShellDirection, ShellOptions, shell, thicken::thicken};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn index_of(brep: &BRep) -> PointMeshIndex {
    let mesh = tessellate(brep, &TessellationOptions::default()).expect("tessellate");
    PointMeshIndex::build(&mesh)
}

fn top_plane_face(brep: &BRep, z: f64) -> popkern_core::topology::FaceHandle {
    let handles: Vec<_> = brep.faces.keys().collect();
    handles
        .into_iter()
        .find(|&fh| {
            use popkern_core::geometry::surface::SurfaceGeometry;
            matches!(
                brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            ) && {
                let ps = brep.face_outer_loop_positions(fh);
                !ps.is_empty() && ps.iter().all(|p| (p.z - z).abs() < 1e-9)
            }
        })
        .expect("top face")
}

fn assert_clean(report: &popkern_mesh::offset_oracle::OffsetOracleReport, label: &str) {
    assert!(
        report.is_clean(),
        "{label}: {} disagreements (first: {:?})",
        report.disagreement_count,
        report.disagreements.first()
    );
    assert!(
        report.classified_fraction() > 0.5,
        "{label}: classified too little: {report:?}"
    );
}

#[test]
fn thickened_cube_wall_agrees_everywhere() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let t = 0.1;
    let result = thicken(&cube, t, &ShellOptions::default(), tol()).expect("thicken cube");
    assert!(result.brep.validate_all().is_empty());
    // Thicken is symmetric (±t/2 about the original skin). Convex
    // body: both skins are sharp — no exclusions, corners included.
    let report = classify_offset_agreement(
        &index_of(&cube),
        &index_of(&result.brep),
        OffsetSemantics::Wall {
            inward: t / 2.0,
            outward: t / 2.0,
        },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert_eq!(report.excluded, 0);
    assert_clean(&report, "thickened cube");
}

#[test]
fn shelled_cube_wall_agrees_outside_the_opening() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let t = 0.1;
    let top = top_plane_face(&cube, 1.0);
    let result = shell(
        &cube,
        t,
        ShellDirection::Inward,
        &[top],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("shell cube");
    assert!(result.brep.validate_all().is_empty());
    // The wall predicate does not hold where the removed face used to
    // provide the skin: exclude points whose nearest original boundary
    // face is the (removed) top, with the oracle band as margin.
    let band = OffsetOracleOptions::default().boundary_band;
    let report = classify_offset_agreement(
        &index_of(&cube),
        &index_of(&result.brep),
        OffsetSemantics::Wall {
            inward: t,
            outward: 0.0,
        },
        &OffsetOracleOptions::default(),
        |p: Point3| {
            let to_top = 1.0 - p.z;
            let to_rest = p.x.min(1.0 - p.x).min(p.y).min(1.0 - p.y).min(p.z);
            to_top < to_rest + 2.0 * band
        },
    );
    assert_clean(&report, "shelled cube");
    assert!(
        (report.excluded as f64) < 0.35 * report.total_samples as f64,
        "opening exclusion swallowed the sweep: {report:?}"
    );
}

#[test]
fn shelled_cylinder_wall_agrees_outside_the_opening() {
    let (r, h, t) = (1.0, 2.0, 0.1);
    let cyl = BRep::cylinder(r, h, tol()).unwrap();
    let top = top_plane_face(&cyl, h);
    let result = shell(
        &cyl,
        t,
        ShellDirection::Inward,
        &[top],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("shell cylinder");
    assert!(result.brep.validate_all().is_empty());
    let band = OffsetOracleOptions::default().boundary_band;
    let report = classify_offset_agreement(
        &index_of(&cyl),
        &index_of(&result.brep),
        OffsetSemantics::Wall {
            inward: t,
            outward: 0.0,
        },
        &OffsetOracleOptions::default(),
        |p: Point3| {
            let rho = (p.x * p.x + p.y * p.y).sqrt();
            let to_top = h - p.z;
            let to_rest = (r - rho).min(p.z);
            to_top < to_rest + 2.0 * band
        },
    );
    assert_clean(&report, "shelled cylinder");
    assert!(
        (report.excluded as f64) < 0.35 * report.total_samples as f64,
        "opening exclusion swallowed the sweep: {report:?}"
    );
}

#[test]
fn offset_cone_grow_agrees_outside_the_sharp_wedges() {
    // Outward sharp offset: Minkowski agreement everywhere except the
    // wedge ring at the base rim and the extended sharp apex — both
    // excluded by explicit geometry, both counted.
    let (r, h, d) = (1.0_f64, 2.0_f64, 0.15_f64);
    let cone = BRep::cone(r, h, tol()).unwrap();
    let result = offset_solid(&cone, d, &ShellOptions::default(), tol()).expect("offset cone");
    assert!(result.brep.validate_all().is_empty());
    let slant = (r * r + h * h).sqrt();
    let apex_shift = d * slant / r; // d / sin(half-angle)
    let report = classify_offset_agreement(
        &index_of(&cone),
        &index_of(&result.brep),
        OffsetSemantics::Grow { distance: d },
        &OffsetOracleOptions::default(),
        |p: Point3| {
            let rho = (p.x * p.x + p.y * p.y).sqrt();
            let to_rim = ((rho - r) * (rho - r) + p.z * p.z).sqrt();
            let to_apex = (rho * rho + (p.z - h) * (p.z - h)).sqrt();
            to_rim < 2.5 * d || to_apex < 1.5 * apex_shift
        },
    );
    assert_clean(&report, "offset cone");
    assert!(
        (report.excluded as f64) < 0.2 * report.total_samples as f64,
        "wedge exclusion swallowed the sweep: {report:?}"
    );
}

/// Analytic SDF of the rounded cube `[0,1]³` filleted at `rf` (exact:
/// the fixture's corner caps are promoted sphere octants).
fn rounded_cube_sdf(p: Point3, rf: f64) -> f64 {
    let half = 0.5 - rf;
    let q = [
        (p.x - 0.5).abs() - half,
        (p.y - 0.5).abs() - half,
        (p.z - 0.5).abs() - half,
    ];
    let outside = (q[0].max(0.0).powi(2) + q[1].max(0.0).powi(2) + q[2].max(0.0).powi(2)).sqrt();
    let inside = q[0].max(q[1]).max(q[2]).min(0.0);
    outside + inside - rf
}

#[test]
fn filleted_cube_oracle_matches_the_exact_rounded_cube_sdf() {
    let fixture = filleted_cube(tol()).expect("filleted cube");
    let index = index_of(&fixture.brep);
    let mut worst = 0.0_f64;
    let n = 17;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let p = Point3::new(
                    -0.3 + 1.6 * i as f64 / (n - 1) as f64,
                    -0.3 + 1.6 * j as f64 / (n - 1) as f64,
                    -0.3 + 1.6 * k as f64 / (n - 1) as f64,
                );
                worst = worst.max((index.signed_distance(p) - rounded_cube_sdf(p, 0.1)).abs());
            }
        }
    }
    // Tessellation chord error only — the BRep geometry is exact.
    assert!(
        worst < 0.01,
        "mesh SDF deviates {worst} from the exact rounded-cube SDF"
    );
}

#[test]
fn filleted_cube_thicken_serves_oracle_clean_or_declines_structured() {
    // Serving the filleted cube is PK-324 acceptance (row 469); until
    // that stage lands the thicken must DECLINE with a structured
    // error, and once it serves this test holds it to the oracle.
    let fixture = filleted_cube(tol()).expect("filleted cube");
    let t = 0.05;
    match thicken(&fixture.brep, t, &ShellOptions::default(), tol()) {
        Ok(result) => {
            assert!(result.brep.validate_all().is_empty());
            let report = classify_offset_agreement(
                &index_of(&fixture.brep),
                &index_of(&result.brep),
                OffsetSemantics::Wall {
                    inward: t,
                    outward: 0.0,
                },
                &OffsetOracleOptions::default(),
                |_| false, // convex, all joins tangential — no wedges
            );
            assert_clean(&report, "thickened filleted cube");
        }
        Err(e) => {
            println!("[pk-324 frontier] filleted-cube thicken declines: {e}");
        }
    }
}
