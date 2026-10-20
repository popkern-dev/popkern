//! "Export this section" (PK-856): the built document's exact section
//! by the page's clipping plane, per body, as SVG and DXF.

use popkern_io::section_drawing::SectionDrawing;
use popkern_live::{HostOptions, Session};

#[test]
fn a_bored_plates_section_is_its_ring_drawn_and_measured() {
    let mut session = Session::new(HostOptions::default());
    let source = r#"let m = model();
let plate = m.profile(xy()).rect(10.0, 10.0).sketch("plate").extrude(blind(2.0));
let bore = m.tool(|t| t.profile(xy()).circle(2.0).sketch("bore").extrude(blind(2.0)));
plate.subtract(bore);
m
"#;
    assert!(session.run(source, "plate.rhai").outcome.is_served());
    // The script sets no unit: document units are SI, the plane z = 1.
    let reply = session.section([0.0, 0.0, 1.0], 1.0).expect("section");
    assert_eq!(reply.bodies.len(), 1);
    let ring = 100.0 - 4.0 * std::f64::consts::PI;
    assert!(
        (reply.bodies[0].area - ring).abs() < 1e-9,
        "{:?}",
        reply.bodies
    );
    let drawing = SectionDrawing::from_svg(&reply.svg).expect("its own SVG");
    assert_eq!(drawing.bodies[0].regions[0].holes.len(), 1, "the bore");
    assert_eq!(reply.dxf.matches("0\nPOLYLINE\n").count(), 2);
    // A plane past the part cuts nothing.
    assert!(
        session
            .section([0.0, 0.0, 1.0], 3.0)
            .expect("section")
            .bodies
            .is_empty()
    );
}

#[test]
fn an_assembly_section_names_each_occurrence_it_cuts() {
    let mut session = Session::new(HostOptions::default());
    let source = r#"let a = assembly("pair");
let block = a.part("block", |m| {
    m.profile(xy()).rect(20.0, 20.0).sketch("s").extrude(blind(10.0));
});
let pin = a.part("pin", |m| {
    m.profile(xy()).circle(2.0).sketch("s").extrude(blind(30.0));
});
a.ground(a.place(block, "block"));
a.place(pin, "pin").translate(-15.0, 0.0, -10.0);
a
"#;
    assert!(session.run(source, "pair.rhai").outcome.is_served());
    let reply = session.section([0.0, 0.0, 1.0], 5.0).expect("section");
    let names: Vec<&str> = reply.bodies.iter().map(|b| b.body.as_str()).collect();
    assert_eq!(names, ["block", "pin"]);
    assert!((reply.bodies[0].area - 400.0).abs() < 1e-9);
    assert!((reply.bodies[1].area - 4.0 * std::f64::consts::PI).abs() < 1e-9);
}

#[test]
fn a_section_before_a_build_declines() {
    let mut session = Session::new(HostOptions::default());
    assert!(session.section([0.0, 0.0, 1.0], 0.0).is_err());
}
