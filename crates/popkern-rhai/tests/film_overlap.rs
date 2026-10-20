//! A film that drives one part through another warns at the frame it
//! does (PK-855); the gallery films, whose fasteners slide in along
//! their own axes and turn where their threads engage, raise none.

use popkern_assembly::Document;
use popkern_rhai::film::{Film, FilmContext, FilmHost, film_name, film_overlaps};
use popkern_rhai::{HostOptions, ScriptHost};

fn document(source: &str, name: &str) -> Document {
    let run = ScriptHost::new()
        .run(source, name)
        .unwrap_or_else(|e| panic!("{name}: {}", e.render()));
    match (run.model, run.assembly) {
        (_, Some(assembly)) => Document::Assembly(Box::new(assembly)),
        (Some(model), None) => Document::Part(Box::new(model)),
        _ => panic!("{name}: no document"),
    }
}

fn film(doc: &mut Document, source: &str, name: &str) -> (Film, f64) {
    let context = FilmContext::from_document(doc);
    let scale = context.unit_scale;
    let run = FilmHost::new(HostOptions::default())
        .run(source, name, &film_name(name), context)
        .unwrap_or_else(|e| panic!("{name}: {}", e.render()));
    (run.film, scale)
}

/// A 20 x 20 x 10 bracket and an r 2 screw standing 3 clear of its
/// -x face.
const BRACKET: &str = r#"let a = assembly("bracket");
let bracket = a.part("bracket", |m| {
    m.profile(xy()).rect(20.0, 20.0).sketch("s").extrude(blind(10.0));
});
let screw = a.part("screw", |m| {
    m.profile(xy()).circle(2.0).sketch("s").extrude(blind(30.0));
});
a.ground(a.place(bracket, "bracket"));
a.place(screw, "screw").translate(-15.0, 0.0, -10.0);
a
"#;

#[test]
fn a_screw_slid_sideways_through_a_bracket_warns_at_the_frame_it_enters() {
    let mut doc = document(BRACKET, "bracket.rhai");
    // 40 to +x over 2 s at 20 fps: 1 a frame. The screw's side reaches
    // the bracket at offset 3, frame 3.
    let (f, scale) = film(
        &mut doc,
        r#"let film = storyboard(#{ fps: 20 });
film.step("Through", "", 2.0, |s| {
    s.move("screw", #{ to: [40.0, 0.0, 0.0], ease: "linear" });
});
film"#,
        "films/through.rhai",
    );
    let Document::Assembly(assembly) = &mut doc else {
        panic!("an assembly")
    };
    let overlaps = film_overlaps(assembly, &f, scale).expect("checked");
    assert_eq!(overlaps.len(), 1, "{overlaps:?}");
    let o = &overlaps[0];
    assert_eq!(o.pair, ["bracket".to_string(), "screw".to_string()]);
    assert!(
        (3..=5).contains(&o.frame),
        "first overlap at frame {}",
        o.frame
    );
    assert_eq!(o.line, Some(3), "the move's line");
    // The same move kept clear of the bracket (below it) raises none.
    let (clear, _) = film(
        &mut doc,
        r#"let film = storyboard(#{ fps: 20 });
film.step("Under", "", 2.0, |s| {
    s.move("screw", #{ to: [0.0, 0.0, -25.0], duration: 1.0 });
    s.move("screw", #{ to: [40.0, 0.0, 0.0], at: 1.0, duration: 1.0, ease: "linear" });
});
film"#,
        "films/under.rhai",
    );
    let Document::Assembly(assembly) = &mut doc else {
        panic!("an assembly")
    };
    assert!(
        film_overlaps(assembly, &clear, scale)
            .expect("checked")
            .is_empty()
    );
}

#[test]
fn the_gallery_assembly_films_drive_no_part_through_another() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");
    let mut films: Vec<std::path::PathBuf> = std::fs::read_dir(&dir)
        .expect("examples")
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.to_str().is_some_and(|s| s.ends_with(".film.rhai")))
        .collect();
    films.sort();
    let mut checked = 0;
    for path in films {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_string();
        let stem = name.split('.').next().unwrap_or_default();
        let film_source = std::fs::read_to_string(&path).expect("film");
        // Only moves and turns (and the insertions recording them) can
        // drive a part into another (a joint
        // sweep is the mates' own swing check): the arm's joint film is
        // skipped without building the arm.
        if ![".move(", ".turn(", ".insert(", ".screw_in("]
            .iter()
            .any(|motion| film_source.contains(motion))
        {
            continue;
        }
        let script = dir.join(format!("{stem}.rhai"));
        let source = std::fs::read_to_string(&script).expect("script");
        let mut doc = document(&source, &format!("{stem}.rhai"));
        if !matches!(doc, Document::Assembly(_)) {
            continue;
        }
        let (f, scale) = film(&mut doc, &film_source, &name);
        let Document::Assembly(assembly) = &mut doc else {
            unreachable!()
        };
        let overlaps = film_overlaps(assembly, &f, scale).expect("checked");
        assert!(overlaps.is_empty(), "{name}: {overlaps:#?}");
        // Control: the nut slid along the bolt WITHOUT its turns drives
        // its thread through the bolt's, and the check sees it.
        if stem == "fastener_kit" {
            let unturned: String = film_source
                .lines()
                .filter(|l| !l.contains("s.turn(\"nut_1\""))
                .collect::<Vec<_>>()
                .join("\n");
            let (f, scale) = film(&mut doc, &unturned, &name);
            let Document::Assembly(assembly) = &mut doc else {
                unreachable!()
            };
            let overlaps = film_overlaps(assembly, &f, scale).expect("checked");
            assert!(
                overlaps.iter().any(|o| o.pair.iter().any(|p| p == "nut_1")),
                "the unturned nut must warn: {overlaps:?}"
            );
        }
        checked += 1;
    }
    // The angled fitting, the fastener kit and the planetary gearbox
    // move or turn parts; the others are part documents or sweep joints.
    assert!(checked >= 3, "{checked} assembly films checked");
}
