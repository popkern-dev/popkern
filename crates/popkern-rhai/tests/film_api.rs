//! Presentation films' scripting surface (docs/architecture/
//! presentation-films.md §10): the storyboard's data, SI lengths from
//! a mm document, and every validation error at its own line.

use popkern_assembly::Document;
use popkern_rhai::film::{Channel, FilmContext, FilmHost, SectionKey, film_name};
use popkern_rhai::{HostOptions, ScriptHost, ScriptRunError};

fn document(source: &str) -> Document {
    let run = ScriptHost::new()
        .run(source, "doc.rhai")
        .expect("document runs");
    match (run.model, run.assembly) {
        (Some(model), _) => Document::Part(Box::new(model)),
        (None, Some(assembly)) => Document::Assembly(Box::new(assembly)),
        _ => panic!("no document"),
    }
}

const SLIDE: &str = r#"let a = assembly("slide");
a.param("gap", 0.5, "m");
a.param("width", 1.0, "m");
let block = a.part("block", |m| {
    m.param("width", 1.0, "m");
    m.profile(xy()).rect(p("width"), 1.0).sketch("s").extrude(blind(1.0));
});
a.drive("width", [block]);
let base = a.place(block, "base");
let top = a.place(block, "top");
a.ground(base);
a.distance(top.datum("xy"), base.datum("xy"), p("gap"), #{ aligned: true });
a.concentric(top.datum("z"), base.datum("z"));
let pair = a.subassembly("pair", |s| { s.place(block, "low"); s.place(block, "high").translate(0.0, 0.0, 2.0); });
a.place(pair).translate(3.0, 0.0, 0.0);
a
"#;

fn film(doc: &mut Document, source: &str) -> Result<popkern_rhai::film::Film, ScriptRunError> {
    let context = FilmContext::from_document(doc);
    FilmHost::new(HostOptions::default())
        .run(source, "guide.rhai", "guide", context)
        .map(|run| run.film)
}

#[test]
fn an_assembly_film_is_steps_and_tracks_at_absolute_times() {
    let mut doc = document(SLIDE);
    let f = film(
        &mut doc,
        r#"let doc = document();
let film = storyboard(#{ size: [1280, 720], fps: 24, look: "studio" });
film.title("Slide", "two blocks and a pair", 2.0, |s| { s.orbit(30.0); });
film.step("Lift", "the top block rises", 3.0, |s| {
    s.move("top", #{ from: [0.0, 0.0, 2.0], at: 0.5, duration: 2.0, ease: "out" });
    s.joint("gap", 0.5, 1.5);
    s.camera(#{ zoom: 1.5, target: "top" });
});
film.hold(1.0);
film.step("Inside", 2.0, |s| {
    s.section(#{ from: "x", from_at: 1.0, to_at: 0.4 });
    s.fade("pair_1", 1.0, 0.3);
    s.show("pair_1/high", false, #{ at: 1.0 });
    s.explode(0.0, 0.5);
});
print(doc.instances().len());
film
"#,
    )
    .expect("film runs");
    assert_eq!(f.name, "guide");
    assert_eq!(f.size, [1280, 720]);
    assert_eq!(f.fps, 24);
    assert_eq!(f.look, "studio");
    assert_eq!(f.duration, 8.0);
    let kinds: Vec<_> = f
        .steps
        .iter()
        .map(|s| (s.kind, s.start, s.duration))
        .collect();
    assert_eq!(
        kinds,
        [
            ("title", 0.0, 2.0),
            ("step", 2.0, 3.0),
            ("hold", 5.0, 1.0),
            ("step", 6.0, 2.0)
        ]
    );
    assert_eq!(f.steps[1].line, Some(4));
    // Times are absolute; a move's at/duration lands inside its step.
    let mv = f
        .tracks
        .iter()
        .find(|t| matches!(t.channel, Channel::Move { .. }))
        .expect("move");
    assert_eq!((mv.start, mv.duration, mv.line), (2.5, 2.0, Some(5)));
    let Channel::Move { target, from, to } = &mv.channel else {
        unreachable!()
    };
    assert_eq!(
        (target.as_str(), *from, *to),
        ("top", [0.0, 0.0, 2.0], [0.0; 3])
    );
    // An orbit is linear by default; others ease in and out.
    assert_eq!(f.tracks[0].ease, popkern_rhai::film::Ease::Linear);
    let section = f
        .tracks
        .iter()
        .find(|t| matches!(t.channel, Channel::Section { .. }))
        .expect("section");
    let Channel::Section { from, to, .. } = &section.channel else {
        unreachable!()
    };
    assert_eq!(
        from,
        &SectionKey::Axis {
            axis: "x".into(),
            at: 1.0
        }
    );
    assert_eq!(
        to,
        &SectionKey::Axis {
            axis: "x".into(),
            at: 0.4
        }
    );
    let show = f
        .tracks
        .iter()
        .find(|t| matches!(t.channel, Channel::Show { .. }))
        .expect("show");
    assert_eq!((show.start, show.duration), (7.0, 0.0));
    // The JSON the page reads.
    let json = serde_json::to_value(&f).expect("json");
    assert_eq!(json["tracks"][0]["channel"]["kind"], "orbit");
    assert_eq!(json["camera"]["target"], "all");
    let camera = &f
        .tracks
        .iter()
        .find(|t| matches!(t.channel, Channel::Camera { .. }))
        .expect("camera");
    let cj = serde_json::to_value(&camera.channel).expect("json");
    assert!(
        cj["to"]["azimuth"].is_null() && cj["to"]["zoom"] == 1.5 && cj["to"]["target"] == "top",
        "{cj}"
    );
}

#[test]
fn a_mm_part_film_moves_in_si_metres() {
    let mut doc = document(
        r#"let m = model(#{ unit: "mm" });
m.profile(xy()).rect(10.0, 10.0).sketch("a").extrude(blind(5.0));
m.profile(xy_at(20.0)).rect(10.0, 10.0).sketch("b").extrude(blind(5.0));
m
"#,
    );
    let f = film(
        &mut doc,
        r#"let film = storyboard();
film.step("Up", 1.0, |s| {
    s.move("body 1", #{ to: [0.0, 0.0, 10.0] });
    s.turn("body 0", axis_at([5.0, 0.0, 0.0], [0.0, 0.0, 2.0]), 0.0, 90.0);
    s.section(#{ from: #{ normal: [0.0, 0.0, 1.0], offset: 2.5 } });
});
film
"#,
    )
    .expect("film runs");
    let Channel::Move { to, .. } = &f.tracks[0].channel else {
        panic!("move")
    };
    assert!((to[2] - 0.01).abs() < 1e-15, "10 mm is 0.01 m: {to:?}");
    let Channel::Turn {
        origin, direction, ..
    } = &f.tracks[1].channel
    else {
        panic!("turn")
    };
    assert!(
        (origin[0] - 0.005).abs() < 1e-15 && direction == &[0.0, 0.0, 1.0],
        "{origin:?} {direction:?}"
    );
    let Channel::Section { from, .. } = &f.tracks[2].channel else {
        panic!("section")
    };
    assert_eq!(
        from,
        &SectionKey::Plane {
            normal: [0.0, 0.0, 1.0],
            offset: 0.0025
        }
    );
}

fn error(doc: &mut Document, source: &str) -> ScriptRunError {
    film(doc, source).expect_err("the film declines")
}

#[test]
fn unknown_names_fail_at_their_line_listing_what_exists() {
    let mut doc = document(SLIDE);
    let e = error(
        &mut doc,
        "let film = storyboard();\nfilm.step(\"a\", 1.0, |s| {\n    s.move(\"lid\", #{ from: [0.0, 0.0, 1.0] });\n});\nfilm\n",
    );
    assert_eq!(e.kind, "ArgumentError", "{e}");
    assert_eq!(e.line, Some(3), "{e}");
    assert!(
        e.message.contains("no 'lid'") && e.message.contains("pair_1/high"),
        "{e}"
    );
    // A prefix selects what lies under it; a partial name does not.
    assert!(film(&mut doc, "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.fade(\"pair_1\", 1.0, 0.0); }); film").is_ok());
    let e = error(
        &mut doc,
        "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.fade(\"pair\", 1.0, 0.0); }); film",
    );
    assert!(e.message.contains("no 'pair'"), "{e}");
}

#[test]
fn a_joint_is_a_placement_only_assembly_parameter() {
    let mut doc = document(SLIDE);
    let e = error(
        &mut doc,
        "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.joint(\"width\", 1.0, 2.0); }); film",
    );
    assert!(e.message.contains("drives a part's geometry"), "{e}");
    let e = error(
        &mut doc,
        "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.joint(\"swing\", 0.0, 1.0); }); film",
    );
    assert!(e.message.contains("joints are gap"), "{e}");
    let mut part = document(
        "let m = model(); m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(1.0)); m",
    );
    let e = error(
        &mut part,
        "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.joint(\"gap\", 0.0, 1.0); }); film",
    );
    assert!(e.message.contains("part document has none"), "{e}");
    let e = error(
        &mut part,
        "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.explode(0.0, 1.0); }); film",
    );
    assert!(e.message.contains("part document has none"), "{e}");
    let e = error(
        &mut part,
        "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.move(\"body 3\", #{ to: [1.0, 0.0, 0.0] }); }); film",
    );
    assert!(
        e.message.contains("no 'body 3'") && e.message.contains("body 0"),
        "{e}"
    );
}

#[test]
fn malformed_storyboards_decline_by_name() {
    let mut doc = document(SLIDE);
    let cases: &[(&str, &str)] = &[
        ("let x = 1; x", "returns its storyboard"),
        ("let film = storyboard(); film", "no steps"),
        (
            "let film = storyboard(); film.step(\"a\", 0.0, |s| {}); film",
            "longer than 0",
        ),
        (
            "let film = storyboard(); film.step(\"a\", 900.0, |s| {}); film",
            "seconds between",
        ),
        (
            "let film = storyboard(#{ fps: 500 }); film.hold(1.0); film",
            "1..120",
        ),
        (
            "let film = storyboard(#{ size: [8000, 100] }); film.hold(1.0); film",
            "sides between",
        ),
        (
            "let film = storyboard(#{ look: \"neon\" }); film.hold(1.0); film",
            "\"toon\"",
        ),
        (
            "let film = storyboard(#{ background: \"white\" }); film.hold(1.0); film",
            "#rrggbb",
        ),
        (
            "let film = storyboard(#{ fsp: 30 }); film.hold(1.0); film",
            "no option 'fsp'",
        ),
        (
            "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.orbit(10.0, #{ ease: \"bounce\" }); }); film",
            "\"linear\"",
        ),
        (
            "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.move(\"top\", #{}); }); film",
            "from: [x, y, z]",
        ),
        (
            "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.fade(\"top\", 1.0, 2.0); }); film",
            "between 0 and 1",
        ),
        (
            "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.section(#{ from: \"w\" }); }); film",
            "\"x\", \"y\", \"z\"",
        ),
        (
            "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.camera(#{ elevation: 95.0 }); }); film",
            "-89..89",
        ),
        (
            "let film = storyboard(); film.step(\"a\", 1.0, |s| { s.turn(\"top\", axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]), 0.0, 90.0); }); film",
            "non-zero",
        ),
    ];
    for (source, needle) in cases {
        let e = error(&mut doc, source);
        assert!(e.to_string().contains(needle), "{source}\n  -> {e}");
    }
}

#[test]
fn a_runaway_film_script_hits_the_operation_budget() {
    let mut doc = document(SLIDE);
    let context = FilmContext::from_document(&mut doc);
    let options = HostOptions {
        max_operations: Some(10_000),
        ..HostOptions::default()
    };
    let e = FilmHost::new(options)
        .run(
            "let film = storyboard(); loop { } film",
            "loop.rhai",
            "loop",
            context,
        )
        .expect_err("budget");
    assert_eq!(e.kind, "BudgetExceeded", "{e}");
}

#[test]
fn a_film_reads_the_documents_parameter_values() {
    let mut doc = document(SLIDE);
    let f = film(
        &mut doc,
        "let doc = document(); let film = storyboard(); film.step(\"a\", 1.0, |s| { s.joint(\"gap\", doc.value(\"gap\"), doc.value(\"gap\") * 3.0); }); film",
    )
    .expect("film");
    let Channel::Joint { from, to, .. } = &f.tracks[0].channel else {
        panic!("joint")
    };
    assert!(
        (from - 0.5).abs() < 1e-12 && (to - 1.5).abs() < 1e-12,
        "{from} {to}"
    );
    let e = error(
        &mut doc,
        "let doc = document(); doc.value(\"nope\"); let film = storyboard(); film.hold(1.0); film",
    );
    assert!(
        e.message.contains("no parameter 'nope'") && e.message.contains("gap"),
        "{e}"
    );
}

#[test]
fn film_names_come_from_their_files() {
    assert_eq!(film_name("films/assembly_guide.rhai"), "assembly_guide");
    assert_eq!(film_name("angled_fitting.guide.film.rhai"), "guide");
    assert_eq!(film_name("angled_fitting.film.rhai"), "angled_fitting");
    assert_eq!(film_name("turntable"), "turntable");
}
