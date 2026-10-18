//! The ball bearing gallery races: bore lead-ins and groove-rim
//! fillets (PK-433).
//!
//! Four defects stood between the gallery's races and a filleted
//! groove rim, each measured against an analytic or OpenCASCADE
//! volume:
//!
//! - a chamfer on a REVOLVED bore went outward: the annulus cap is
//!   bounded by two one-vertex circles, its vertex Newell is zero, and
//!   the interior hint stepped toward the boundary centroid — which
//!   lies in the hole — so the bore lead-ins added a ring of material
//!   with validators clean;
//! - the groove torus and the wall cylinder BOTH carry a seam at the
//!   rim vertex, and the closed-rim surgery handled one seamed side:
//!   the second rim's torus loop was orphaned (broken half-edge cycle);
//! - a G1 rim fillet between the coaxial torus and cylinder fell to
//!   the generic arc sweep, which exported as a degenerate strip
//!   (OpenCASCADE area 3e-7) — it is exactly a torus;
//! - the G2 strip's seam ends projected to opposite clamped ends of
//!   the closed NURBS chart and the seam ran diagonally through the
//!   antipode (pcurves 1.07e-2 off, OCC: self-intersecting wire).

use popkern_core::topology::BRep;

const RACE: &str = "let m = model();
let pitch = 0.7; let groove_r = 0.16; let lip = sqrt(groove_r * groove_r - 0.12 * 0.12);";

const INNER_SECTION: &str = r#"
let race = m.profile(xz())
    .path(0.35, 0.0).line_to(0.58, 0.0).line_to(0.58, 0.2 - lip)
    .arc_to([pitch - groove_r, 0.2], [0.58, 0.2 + lip])
    .line_to(0.58, 0.4).line_to(0.35, 0.4).close()
    .sketch("s").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());"#;

const OUTER_SECTION: &str = r#"
let race = m.profile(xz())
    .path(0.82, 0.2 - lip).line_to(0.82, 0.0).line_to(1.05, 0.0)
    .line_to(1.05, 0.4).line_to(0.82, 0.4).line_to(0.82, 0.2 + lip)
    .arc_to([pitch + groove_r, 0.2], [0.82, 0.2 - lip]).close()
    .sketch("o").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());"#;

const LEAD_INS: &str = "
race = race.chamfer(edges().circular().at_height(0.0).or(edges().circular().at_height(0.4)), 0.02);";

/// The two rims where the groove torus meets the race wall.
fn groove_rims(continuity: &str) -> String {
    format!(
        "\nrace = race.fillet(edges().circular().at_height(0.2 - lip).or(edges().circular().at_height(0.2 + lip)), 0.012, #{{ continuity: \"{continuity}\" }});"
    )
}

fn build(src: &str) -> BRep {
    let mut run = popkern_rhai::ScriptHost::new()
        .run(&format!("{src}\nm"), "race.rhai")
        .expect("race script serves");
    let parts = run
        .model
        .take()
        .expect("model")
        .to_parts_all()
        .expect("parts");
    let body = parts[0].brep().clone();
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    body
}

fn volume(body: &BRep) -> f64 {
    popkern_measure::body_volume(body, &body.tolerance).expect("volume")
}

fn surface_kinds(body: &BRep) -> Vec<&'static str> {
    body.faces
        .values()
        .filter_map(|f| f.surface.as_ref().map(|s| s.geometry.kind_name()))
        .collect()
}

/// Pappus: the race's revolved section, groove segment and four 0.02
/// lead-in triangles.
fn analytic_chamfered(r_in: f64, r_out: f64, groove_inward: bool) -> f64 {
    use std::f64::consts::PI;
    let annulus = PI * (r_out * r_out - r_in * r_in) * 0.4;
    let (r, half) = (0.16_f64, 0.75_f64.acos());
    let lens = half - half.sin() * half.cos();
    let centroid = 2.0 * r * half.sin().powi(3) / (3.0 * lens);
    let groove_r = if groove_inward {
        0.7 - centroid
    } else {
        0.7 + centroid
    };
    let groove = r * r * lens * 2.0 * PI * groove_r;
    let lead_ins = 2.0 * 0.0002 * 2.0 * PI * ((r_in + 0.02 / 3.0) + (r_out - 0.02 / 3.0));
    annulus - groove - lead_ins
}

#[test]
fn a_revolved_bore_lead_in_cuts_into_the_race() {
    for (section, r_in, r_out, inward) in [
        (INNER_SECTION, 0.35, 0.58, true),
        (OUTER_SECTION, 0.82, 1.05, false),
    ] {
        let body = build(&format!("{RACE}{section}{LEAD_INS}"));
        let want = analytic_chamfered(r_in, r_out, inward);
        let got = volume(&body);
        assert!((got - want).abs() < 1e-9, "volume {got}, analytic {want}");
    }
}

/// Exact Pappus volumes of the G1-filleted races: pi * (closed
/// integral of r^2 dz) round the section, whose fillet arcs are
/// centred where the wall's and the groove's +-0.012 offsets cross.
/// (OpenCASCADE's own fillets of the same races land 7e-9 and 1.6e-8
/// off; its measure of popkern's exported torus body agrees with these
/// to 1e-14.)
const EXACT_G1_INNER: f64 = 0.245_885_947_709_960_58;
const EXACT_G1_OUTER: f64 = 0.505_271_789_848_622_1;

#[test]
fn g1_groove_rims_serve_exact_tori_at_the_pappus_volume() {
    for (section, occ) in [
        (INNER_SECTION, EXACT_G1_INNER),
        (OUTER_SECTION, EXACT_G1_OUTER),
    ] {
        let body = build(&format!("{RACE}{section}{LEAD_INS}{}", groove_rims("g1")));
        let kinds = surface_kinds(&body);
        assert_eq!(kinds.len(), 12, "{kinds:?}");
        assert_eq!(
            kinds.iter().filter(|k| **k == "Torus").count(),
            3,
            "{kinds:?}"
        );
        let got = volume(&body);
        assert!((got - occ).abs() < 1e-10, "volume {got}, exact {occ}");
    }
}

#[test]
fn g2_groove_rims_serve_revolved_strips_fuller_than_the_ball() {
    // A G2 section leaves each support with the support's own section
    // curvature (zero on the wall), so it runs outside the rolling ball
    // and keeps a little more of the corner than G1. The races are
    // mirror images in the meridian plane: the same excess section,
    // revolved at radii in the ratio of the walls. The lofted strip
    // broke that symmetry (its rails sagged off the contact circles
    // between stations: inner excess 3e-9, outer 7.4e-6).
    let excess = |section: &str, exact: f64| {
        let body = build(&format!("{RACE}{section}{LEAD_INS}{}", groove_rims("g2")));
        assert_eq!(body.faces.len(), 12);
        volume(&body) - exact
    };
    let inner = excess(INNER_SECTION, EXACT_G1_INNER);
    let outer = excess(OUTER_SECTION, EXACT_G1_OUTER);
    assert!(inner > 1e-6 && inner < 2e-5, "inner G2 excess {inner}");
    let ratio = outer / inner;
    assert!(
        (ratio - 0.82 / 0.58).abs() < 0.02,
        "outer/inner excess {ratio}"
    );
}
