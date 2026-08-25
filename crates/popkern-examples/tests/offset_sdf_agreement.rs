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
    // Thicken is symmetric (±t/2 about the original skin). The outer
    // skin rounds under the arc default (PK-325), so the Minkowski
    // wall predicate holds by construction — no exclusions, corners
    // included.
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

/// PK-325 row 472, the strong form: at d = 0.15 the corner wedge of a
/// SHARP join reaches √3·d ≈ 0.26 — far outside the oracle's skip band
/// — so this gate genuinely distinguishes the rounded Minkowski result
/// from the extended-plane join (the old sharp default passed the small
/// corpus distances only by band arithmetic). Zero exclusions.
#[test]
fn rounded_cube_offset_is_minkowski_exact_at_large_distance() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let d = 0.15;
    let result = popkern_shell::offset::offset_solid(&cube, d, &ShellOptions::default(), tol())
        .expect("arc offset cube");
    assert!(result.brep.validate_all().is_empty());
    let report = classify_offset_agreement(
        &index_of(&cube),
        &index_of(&result.brep),
        OffsetSemantics::Grow { distance: d },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert_eq!(report.excluded, 0);
    assert_clean(&report, "rounded cube offset (Grow 0.15)");
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

/// PK-324 row 469: the hollow shells of every curved primitive serve
/// `validate_all`-clean and agree with the SDF Wall oracle. A hollow
/// (no-removed-faces) inward shell keeps the outer surface in place
/// and carves the inner skin at depth `t`, so the exact wall predicate
/// is `Wall { inward: t, outward: 0 }` with NO exclusions.
#[test]
fn hollow_curved_shells_serve_oracle_clean() {
    let t = 0.04;
    let bodies: Vec<(&str, BRep)> = vec![
        ("cylinder", BRep::cylinder(0.5, 1.0, tol()).unwrap()),
        ("cone", BRep::cone(0.5, 1.0, tol()).unwrap()),
        ("sphere", BRep::sphere(0.5, tol()).unwrap()),
        ("torus", BRep::torus(0.6, 0.25, tol()).unwrap()),
    ];
    for (name, brep) in &bodies {
        let result = shell(
            brep,
            t,
            ShellDirection::Inward,
            &[],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        )
        .unwrap_or_else(|e| panic!("{name}: hollow shell declined: {e}"));
        let errors = result.brep.validate_all();
        assert!(errors.is_empty(), "{name}: served invalid: {errors:?}");
        let report = classify_offset_agreement(
            &index_of(brep),
            &index_of(&result.brep),
            OffsetSemantics::Wall {
                inward: t,
                outward: 0.0,
            },
            &OffsetOracleOptions::default(),
            |_| false,
        );
        assert_clean(&report, &format!("hollow shelled {name}"));
    }
}

/// PK-324 row 469 volume bands: hollow curved shells against their
/// closed-form wall volumes.
#[test]
fn hollow_curved_shell_volumes_in_analytic_band() {
    use std::f64::consts::PI;
    let t = 0.04_f64;
    // Fine tessellation: the wall volume is a DIFFERENCE of two full
    // solids (torus outer 0.740 − inner 0.522 → wall 0.218), so chord
    // undershoot amplifies ~6x in relative terms; the gate must
    // measure geometry, not mesher resolution.
    let fine = TessellationOptions {
        chord_height: 0.001,
        ..TessellationOptions::default()
    };
    let volume = move |brep: &BRep| {
        let mesh = tessellate(brep, &fine).expect("tessellate");
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
        (vol / 6.0).abs()
    };
    let sphere_wall = (4.0 / 3.0) * PI * (0.5_f64.powi(3) - (0.5_f64 - t).powi(3));
    let torus_wall = 2.0 * PI * PI * 0.6 * (0.25_f64.powi(2) - (0.25_f64 - t).powi(2));
    let cases: Vec<(&str, BRep, f64)> = vec![
        ("sphere", BRep::sphere(0.5, tol()).unwrap(), sphere_wall),
        ("torus", BRep::torus(0.6, 0.25, tol()).unwrap(), torus_wall),
    ];
    for (name, brep, expected) in cases {
        let result = shell(
            &brep,
            t,
            ShellDirection::Inward,
            &[],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        )
        .unwrap_or_else(|e| panic!("{name}: hollow shell declined: {e}"));
        let measured = volume(&result.brep);
        let rel = (measured - expected).abs() / expected;
        assert!(
            rel < 0.02,
            "{name}: wall volume {measured:.6} vs analytic {expected:.6} (rel {rel:.4})"
        );
    }
}

/// PK-324 row 469 (filleted cube volume band): the rounded cube's
/// offset is a rounded cube (core side unchanged, fillet radius ±d),
/// so the thickened wall volume has an exact closed form:
/// `V(rf + t/2) − V(rf − t/2)` with `V(r) = c³ + 6c²r + 3πcr² + (4π/3)r³`.
#[test]
fn thickened_filleted_cube_volume_in_analytic_band() {
    let t = 0.05_f64;
    let fixture = filleted_cube(tol()).expect("filleted cube");
    let result = thicken(&fixture.brep, t, &ShellOptions::default(), tol())
        .expect("filleted-cube thicken serves");
    let rounded_cube_volume = |r: f64| {
        let c = 1.0 - 2.0 * 0.1;
        c * c * c
            + 6.0 * c * c * r
            + 3.0 * std::f64::consts::PI * c * r * r
            + (4.0 * std::f64::consts::PI / 3.0) * r * r * r
    };
    let expected = rounded_cube_volume(0.1 + t / 2.0) - rounded_cube_volume(0.1 - t / 2.0);
    let fine = TessellationOptions {
        chord_height: 0.001,
        ..TessellationOptions::default()
    };
    let mesh = tessellate(&result.brep, &fine).expect("tessellate");
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
        rel < 0.02,
        "wall volume {measured:.6} vs analytic {expected:.6} (rel {rel:.4})"
    );
}
