//! `s.evolve` (PK-959): a generative design converging. The design is
//! named by its feature label; its layer is the selector
//! `"design <label>"`, which show and fade reach like a body.

use popkern_assembly::Document;
use popkern_rhai::film::{Channel, Ease, Film, FilmContext, FilmDesign, FilmHost};
use popkern_rhai::{HostOptions, ScriptHost, ScriptRunError};

/// A part's film context with one design, `carve`, solved in 120
/// updates (the context a generative part's build gives).
fn context() -> FilmContext {
    let run = ScriptHost::new()
        .run(
            "let m = model(); m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(1.0)); m",
            "doc.rhai",
        )
        .expect("document runs");
    let mut doc = Document::Part(Box::new(run.model.expect("a part")));
    let mut context = FilmContext::from_document(&mut doc);
    assert!(context.designs.is_empty(), "the block has no design");
    context.designs.push(FilmDesign {
        label: "carve".into(),
        design: "0".repeat(32),
        updates: 120,
    });
    context
}

fn film(context: FilmContext, steps: &str) -> Result<Film, ScriptRunError> {
    let source = format!("let film = storyboard();\n{steps}\nfilm\n");
    FilmHost::new(HostOptions::default())
        .run(&source, "carve.rhai", "carve", context)
        .map(|run| run.film)
}

#[test]
fn evolve_runs_the_design_over_its_step_linearly_by_default() {
    let f = film(
        context(),
        r#"film.step("Carve", 6.0, |s| { s.evolve("carve"); });
film.step("Settle", 2.0, |s| { s.evolve("design carve", #{ from: 60, to: "end", duration: 1.0, ease: "out" }); s.fade("design carve", 1.0, 0.0); s.show("body 0", true); });"#,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let evolves: Vec<_> = f
        .tracks
        .iter()
        .filter_map(|t| match &t.channel {
            Channel::Evolve { target, from, to } => {
                Some((t.start, t.duration, t.ease, target.as_str(), *from, *to))
            }
            _ => None,
        })
        .collect();
    assert_eq!(
        evolves,
        vec![
            (0.0, 6.0, Ease::Linear, "design carve", 0.0, 120.0),
            (6.0, 1.0, Ease::Out, "design carve", 60.0, 120.0),
        ]
    );
    assert!(
        f.tracks.iter().any(
            |t| matches!(&t.channel, Channel::Fade { target, .. } if target == "design carve")
        )
    );
}

#[test]
fn evolve_declines_unknown_designs_and_updates_at_their_line() {
    let cases = [
        (
            r#"s.evolve("shell");"#,
            "has no generative design 'shell'; it has carve",
        ),
        (r#"s.evolve("carve", #{ to: 121 });"#, "between 0 and 120"),
        (
            r#"s.evolve("carve", #{ from: "start" });"#,
            "an update number or \"end\"",
        ),
        (r#"s.evolve("carve", #{ speed: 2 });"#, "speed"),
        (r#"s.fade("design shell", 1.0, 0.0);"#, "design carve"),
    ];
    for (action, needle) in cases {
        let e = film(
            context(),
            &format!("film.step(\"a\", 1.0, |s| {{\n    {action}\n}});"),
        )
        .expect_err("the film declines");
        assert_eq!(e.kind, "ArgumentError", "{action}: {e}");
        assert_eq!(e.line, Some(3), "{action}: {e}");
        assert!(e.to_string().contains(needle), "{action}\n  -> {e}");
    }
    let mut none = context();
    none.designs.clear();
    let e =
        film(none, "film.step(\"a\", 1.0, |s| { s.evolve(\"carve\"); });").expect_err("no design");
    assert!(
        e.to_string().contains("no generative design to evolve"),
        "{e}"
    );
}
