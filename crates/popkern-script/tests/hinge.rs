//! PK-699: hinges as history features — one feature per component
//! sharing one parameter set owned by the first, parameter edits that
//! replay every component, `.pkn` round trips, `only:`, and living
//! hinges.

use popkern_history::feature_hinge::{HingeComponentDesc, HingeKindDesc, LivingHingeStyleDesc};
use popkern_script::body_hinge::{HingeOptions, LivingHingeOptions};
use popkern_script::model::Model;
use popkern_script::scalar::Scalar;

fn faces(model: &mut Model, body: popkern_script::body::Body) -> usize {
    model.brep_of(body).expect("brep").faces.len()
}

/// Every vertex position, in arena order: a bit-exact fingerprint.
fn vertices(model: &mut Model, body: popkern_script::body::Body) -> Vec<[f64; 3]> {
    let brep = model.brep_of(body).expect("brep");
    brep.vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect()
}

/// The largest |x| and the largest z over the body's vertices.
fn extent(model: &mut Model, body: popkern_script::body::Body) -> (f64, f64) {
    vertices(model, body)
        .iter()
        .fold((0.0_f64, 0.0_f64), |(r, z), p| {
            (r.max(p[0].abs()), z.max(p[2]))
        })
}

#[test]
fn a_butt_hinge_records_three_bodies_sharing_one_parameter_set() {
    let mut m = Model::new();
    let h = m
        .hinge(HingeKindDesc::Butt, 60.0, &HingeOptions::default(), "")
        .expect("hinge");
    let (a, b, pin) = (
        h.leaf_a.expect("a"),
        h.leaf_b.expect("b"),
        h.pin.expect("pin"),
    );
    assert_eq!(m.bodies().len(), 3);
    let labels: Vec<String> = m.features().iter().map(|f| f.label.clone()).collect();
    assert_eq!(labels.len(), 3);
    assert!(labels[1] == format!("{}_b", labels[0]), "{labels:?}");
    assert!(labels[2] == format!("{}_pin", labels[0]), "{labels:?}");
    // One slider set, named after the first component; unset dimensions
    // follow the length (expression-driven), counts and the angle do not.
    let params = m.params();
    let owner = &labels[0];
    let get = |role: &str| {
        params
            .iter()
            .find(|p| p.name == format!("{owner}.{role}"))
            .unwrap_or_else(|| {
                panic!(
                    "{owner}.{role} in {:?}",
                    params.iter().map(|p| &p.name).collect::<Vec<_>>()
                )
            })
    };
    assert!(get("length").expression.is_none());
    assert!(get("knuckles").expression.is_none());
    assert!(get("angle").expression.is_none());
    assert!(get("leaf_width").expression.is_some(), "follows the length");
    assert!(
        params
            .iter()
            .all(|p| !p.name.starts_with(&format!("{}.", labels[1])))
    );
    assert!(faces(&mut m, a) > 0 && faces(&mut m, b) > 0 && faces(&mut m, pin) > 0);

    // A knuckle-count edit replays both leaves (more knuckles, more faces).
    let before = (faces(&mut m, a), faces(&mut m, b));
    m.set_param(&format!("{owner}.knuckles"), 7.0)
        .expect("knuckles");
    let after = (faces(&mut m, a), faces(&mut m, b));
    assert!(
        after.0 > before.0 && after.1 > before.1,
        "{before:?} -> {after:?}"
    );
    // A length edit rescales everything that follows it: the pin's
    // radius (L/30) and length double with L.
    let (r, z) = extent(&mut m, pin);
    m.set_param(&format!("{owner}.length"), 120.0)
        .expect("length");
    let (r2, z2) = extent(&mut m, pin);
    assert!(
        (r2 / r - 2.0).abs() < 1e-9 && (z2 / z - 2.0).abs() < 1e-9,
        "{r},{z} -> {r2},{z2}"
    );
}

#[test]
fn hinges_round_trip_through_pkn() {
    let mut m = Model::new();
    let options = HingeOptions {
        knuckles: Some(Scalar::Literal(7.0)),
        clearance: Some(Scalar::Literal(0.3)),
        ..HingeOptions::default()
    };
    m.hinge(HingeKindDesc::PrintInPlace, 60.0, &options, "lid")
        .expect("hinge");
    m.living_hinge(
        LivingHingeStyleDesc::Slotted,
        10.0,
        40.0,
        &LivingHingeOptions::default(),
        "flex",
    )
    .expect("living");
    let before: Vec<_> = m
        .bodies()
        .into_iter()
        .map(|b| vertices(&mut m, b))
        .collect();
    let bytes = m.to_pkn_bytes().expect("save");
    let mut loaded = Model::from_pkn_bytes(&bytes).expect("load");
    loaded.replay().expect("replay");
    let after: Vec<_> = loaded
        .bodies()
        .into_iter()
        .map(|b| vertices(&mut loaded, b))
        .collect();
    // Arena order may differ after a load; the geometry may not: each
    // body's sorted vertex set, bodies in a canonical order.
    let canon = |mut bodies: Vec<Vec<[f64; 3]>>| {
        for b in &mut bodies {
            b.sort_by(|p, q| p.partial_cmp(q).expect("finite"));
        }
        bodies.sort_by(|p, q| p.partial_cmp(q).expect("finite"));
        bodies
    };
    assert_eq!(before.len(), 3, "two leaves and the living hinge");
    let (before, after) = (canon(before), canon(after));
    assert_eq!(before.len(), after.len());
    for (i, (x, y)) in before.iter().zip(&after).enumerate() {
        assert_eq!(x.len(), y.len(), "body {i} vertex count");
        assert_eq!(x, y, "body {i} bit-identical after a round trip");
    }
    loaded
        .set_param("lid.clearance", 0.4)
        .expect("edit after load");
}

#[test]
fn only_records_one_component_which_owns_the_parameters() {
    let mut m = Model::new();
    let options = HingeOptions {
        only: Some(HingeComponentDesc::LeafB),
        ..HingeOptions::default()
    };
    let h = m
        .hinge(HingeKindDesc::Butt, 60.0, &options, "door")
        .expect("leaf b");
    assert!(h.leaf_a.is_none() && h.pin.is_none());
    assert!(h.leaf_b.is_some());
    assert_eq!(m.bodies().len(), 1);
    assert!(m.params().iter().any(|p| p.name == "door.knuckles"));
    m.set_param("door.angle", 90.0).expect("angle");
    let pin_of_print = Model::new().hinge(
        HingeKindDesc::PrintInPlace,
        60.0,
        &HingeOptions {
            only: Some(HingeComponentDesc::Pin),
            ..HingeOptions::default()
        },
        "",
    );
    assert!(pin_of_print.is_err());
}

#[test]
fn a_failing_hinge_reports_its_typed_decline() {
    let mut m = Model::new();
    let options = HingeOptions {
        knuckles: Some(Scalar::Literal(41.0)),
        gap: Some(Scalar::Literal(1.5)),
        ..HingeOptions::default()
    };
    let err = m
        .hinge(HingeKindDesc::Butt, 60.0, &options, "")
        .expect_err("declines");
    let text = err.to_string();
    assert!(text.contains("41 knuckles"), "{text}");
}

#[test]
fn living_hinges_record_one_body_with_style_parameters_only() {
    let mut m = Model::new();
    m.living_hinge(
        LivingHingeStyleDesc::Straight,
        10.0,
        40.0,
        &LivingHingeOptions::default(),
        "flex",
    )
    .expect("straight");
    let names: Vec<String> = m.params().into_iter().map(|p| p.name).collect();
    assert!(names.contains(&"flex.web".to_string()));
    assert!(!names.iter().any(|n| n.contains("slot")), "{names:?}");
    m.set_param("flex.span", 14.0).expect("span edit replays");
}

/// PK-710: a loaded document keeps its bodies (the body manager is not
/// serialized; load rebuilds it from the feature nodes).
#[test]
fn loaded_documents_keep_their_bodies() {
    use popkern_history::feature_fastener::FastenerKind;
    use popkern_script::body_fastener::{FastenerOptions, ThreadSpecArg};
    let mut m = Model::new();
    m.fastener(
        FastenerKind::ThreadedRod,
        ThreadSpecArg::new(6.0, 1.0),
        10.0,
        FastenerOptions::default(),
        "rod",
    )
    .expect("rod");
    m.hinge(HingeKindDesc::Butt, 60.0, &HingeOptions::default(), "h")
        .expect("hinge");
    let ids = |m: &Model| {
        m.bodies()
            .iter()
            .map(|b| format!("{b:?}"))
            .collect::<Vec<_>>()
    };
    let before = ids(&m);
    assert_eq!(before.len(), 4);
    let mut loaded = Model::from_pkn_bytes(&m.to_pkn_bytes().expect("save")).expect("load");
    assert_eq!(ids(&loaded), before, "bodies come back with the document");
    loaded.replay().expect("replay");
    assert_eq!(ids(&loaded), before, "and survive a replay");
    // A body recorded after the load gets a fresh id.
    loaded
        .hinge(HingeKindDesc::Butt, 30.0, &HingeOptions::default(), "h2")
        .expect("second hinge");
    assert_eq!(loaded.bodies().len(), 7);
}
