//! PK-948: the `simulate` request on a part built in popkern live.
//!
//! These are kernel-suite tests (a few solves on a small bar, seconds):
//! they check the live path end to end. Faces are picked by their
//! stable names, the bar's deflection and heat flow match closed forms,
//! a setup survives a parametric edit (gate PE1), and a missing face is
//! a typed decline. Solver validation is the physics suite's job
//! (docs/architecture/analysis-integration.md).

use popkern_live::protocol::Protocol;
use popkern_live::{HostOptions, Session};
use serde_json::{Value, json};

/// A 40 x 4 x 4 mm bar along x, its length a parameter, in a mm
/// document (lengths stored in metres, PK-462).
const BAR: &str = r#"let m = model(#{ unit: "mm" });
m.param("length", 40.0, "mm");
m.profile(xy()).rect(p("length"), 4.0).sketch("bar_sketch").extrude(blind(4.0));
m
"#;

/// The same bar in a default-unit document (millimetre-magnitude
/// numbers, as the gallery's parts are written).
const BAR_DEFAULT_UNITS: &str = r#"let m = model();
m.profile(xy()).rect(40.0, 4.0).sketch("bar_sketch").extrude(blind(4.0));
m
"#;

/// The stable names of the bar's two end faces (least and greatest x),
/// read off the display mesh as the panel's face picks read them.
fn end_faces(session: &mut Session, script: &str) -> (String, String) {
    let envelope = session.run(script, "bar.rhai");
    assert!(envelope.outcome.is_served(), "{:?}", envelope.outcome);
    let body = &envelope.bodies[0];
    let mut mean_x = vec![(0.0f64, 0usize); body.face_count];
    for (t, &face) in body.face_ids.iter().enumerate() {
        for k in 0..3 {
            let v = body.indices[3 * t + k] as usize;
            mean_x[face as usize].0 += f64::from(body.positions[3 * v]);
            mean_x[face as usize].1 += 1;
        }
    }
    let x = |f: usize| mean_x[f].0 / mean_x[f].1.max(1) as f64;
    let faces: Vec<usize> = (0..body.face_count).collect();
    let low = faces.iter().copied().min_by(|&a, &b| x(a).total_cmp(&x(b)));
    let high = faces.iter().copied().max_by(|&a, &b| x(a).total_cmp(&x(b)));
    let label = |f: Option<usize>| body.face_labels[f.expect("face")].clone();
    let (low, high) = (label(low), label(high));
    assert!(low != "unnamed" && high != "unnamed", "{low} {high}");
    (low, high)
}

fn cantilever(held: &str, loaded: &str) -> String {
    format!(
        r#"units = "mm"
material = "PLA (FDM)"

[faces]
held = {{ name = "{held}" }}
tip = {{ name = "{loaded}" }}

[[support]]
faces = ["held"]
kind = "fixed"

[[load]]
faces = ["tip"]
force_n = 10
direction = [0, 0, -1]
"#
    )
}

fn youngs(material: &str) -> f64 {
    popkern_live::session_simulate::materials()
        .as_array()
        .and_then(|m| m.iter().find(|m| m["name"] == material))
        .and_then(|m| m["youngs_modulus_mpa"].as_f64())
        .expect("material")
}

/// Euler-Bernoulli tip deflection of the 4 x 4 bar of length `l` under
/// 10 N, mm.
fn beam_tip(l: f64) -> f64 {
    let i = 4.0 * 4.0f64.powi(3) / 12.0;
    10.0 * l.powi(3) / (3.0 * youngs("PLA (FDM)") * i)
}

fn solved(result: &popkern_live::session_simulate::SimulateResult) -> &Value {
    assert_eq!(result.summary["status"], "solved", "{:#}", result.summary);
    &result.summary["result"]
}

#[test]
fn a_cantilever_bar_bends_as_the_beam_formula_and_survives_a_length_edit() {
    let mut session = Session::new(HostOptions::default());
    let (held, loaded) = end_faces(&mut session, BAR);
    let setup = cantilever(&held, &loaded);
    let first = session.simulate(&setup).expect("simulate");
    let result = solved(&first);
    let tip = result["max_displacement_mm"].as_f64().expect("deflection");
    let want = beam_tip(40.0);
    // Slender (L/h = 10): shear and the clamp's local flexibility add a
    // few percent to the Euler-Bernoulli value.
    assert!(
        tip > want && tip < 1.08 * want,
        "tip {tip} mm vs beam {want} mm"
    );
    // Every field is per vertex of the mesh it came with.
    assert_eq!(first.meshes.len(), 1);
    let mesh = &first.meshes[0];
    let vertices = mesh.positions.len() / 3;
    assert_eq!(mesh.displacement.len(), 3 * vertices);
    for name in [
        "von_mises_mpa",
        "max_principal_mpa",
        "displacement_mm",
        "zz_error_mpa",
    ] {
        assert_eq!(mesh.fields[name].len(), vertices, "{name}");
    }
    assert!(mesh.indices.iter().all(|&i| (i as usize) < vertices));
    // PE1: a shorter bar, the same setup by name.
    let edited = session.set_param("length", 30.0);
    assert!(edited.outcome.is_served(), "{:?}", edited.outcome);
    let second = session.simulate(&setup).expect("simulate");
    let short = solved(&second)["max_displacement_mm"]
        .as_f64()
        .expect("deflection");
    let ratio = short / tip;
    let cube = (30.0f64 / 40.0).powi(3);
    assert!(
        (ratio - cube).abs() < 0.05 * cube,
        "30 mm / 40 mm deflection {ratio} vs (3/4)^3 = {cube}"
    );
}

#[test]
fn heat_through_an_insulated_bar_rises_as_the_conduction_formula() {
    let mut session = Session::new(HostOptions::default());
    let (cold, hot) = end_faces(&mut session, BAR_DEFAULT_UNITS);
    let setup = format!(
        r#"analysis = "thermal"
units = "mm"
material = "Aluminium 6061-T6"

[faces]
cold = {{ name = "{cold}" }}
hot = {{ name = "{hot}" }}

[[fixed_temperature]]
faces = ["cold"]
temperature_c = 20

[[heat_flux]]
faces = ["hot"]
power_w = 1
"#
    );
    let out = session.simulate(&setup).expect("simulate");
    let result = solved(&out);
    let hot_c = result["max_temperature_c"].as_f64().expect("temperature");
    // A linear field along the bar: dT = P L / (k A), with k from the
    // solver's own material (W/(m K) -> W/(mm K)).
    let k = popkern_analysis::material::preset("Aluminium 6061-T6")
        .ok()
        .and_then(|m| m.conductivity.value())
        .expect("conductivity")
        / 1000.0;
    let want = 20.0 + 40.0 / (k * 16.0);
    assert!(
        (hot_c - want).abs() < 0.01 * (want - 20.0),
        "hot end {hot_c} C vs {want} C"
    );
    assert_eq!(
        out.meshes[0].fields["temperature_c"].len(),
        out.meshes[0].positions.len() / 3
    );
}

#[test]
fn a_face_name_the_part_does_not_have_is_a_typed_decline() {
    let mut session = Session::new(HostOptions::default());
    let (held, _) = end_faces(&mut session, BAR);
    let out = session
        .simulate(&cantilever(&held, "bar_sketch.no_such_face"))
        .expect("simulate");
    assert_eq!(out.summary["status"], "declined", "{:#}", out.summary);
    assert_eq!(out.summary["error"]["kind"], "SelectorUnresolved");
    assert!(
        out.summary["error"]["message"]
            .as_str()
            .is_some_and(|m| m.contains("tip")),
        "{:#}",
        out.summary
    );
    assert!(out.meshes.is_empty());
}

#[test]
fn simulate_and_materials_answer_over_the_protocol_with_field_buffers() {
    let mut protocol = Protocol::new();
    let init = protocol.handle(&json!({ "id": 1, "type": "init", "options": {} }));
    assert_eq!(init.json["ok"], true, "{:#}", init.json);
    let run =
        protocol.handle(&json!({ "id": 2, "type": "run", "source": BAR, "name": "bar.rhai" }));
    assert_eq!(run.json["ok"], true, "{:#}", run.json);
    let materials = protocol.handle(&json!({ "id": 3, "type": "materials" }));
    assert!(
        materials.json["materials"]
            .as_array()
            .is_some_and(|m| m.iter().any(|m| m["name"] == "PLA (FDM)"))
    );
    let mut session = Session::new(HostOptions::default());
    let (held, loaded) = end_faces(&mut session, BAR);
    let reply = protocol
        .handle(&json!({ "id": 4, "type": "simulate", "setup": cantilever(&held, &loaded) }));
    assert_eq!(
        reply.json["simulation"]["status"], "solved",
        "{:#}",
        reply.json
    );
    let paths: Vec<&str> = reply.attachments.iter().map(|a| a.path.as_str()).collect();
    for path in [
        "meshes.0.positions",
        "meshes.0.indices",
        "meshes.0.displacement",
        "meshes.0.fields.von_mises_mpa",
    ] {
        assert!(paths.contains(&path), "{path} missing from {paths:?}");
    }
}
