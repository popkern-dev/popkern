//! Film fastener insertion (PK-854): `s.insert` and `s.screw_in` read a
//! fastener's axis and lead from the built assembly. The angled fitting
//! guide written with them moves every screw exactly as its explicit
//! moves and turns did; an occurrence with no axis declines naming it.

use popkern_assembly::Document;
use popkern_rhai::film::{Film, FilmContext, FilmHost, film_name, occurrence_motion};
use popkern_rhai::{HostOptions, ScriptHost};

fn angled_fitting() -> Document {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");
    let source = std::fs::read_to_string(dir.join("angled_fitting.rhai")).expect("script");
    let run = ScriptHost::new()
        .run(&source, "angled_fitting.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    Document::Assembly(Box::new(run.assembly.expect("an assembly")))
}

fn film(doc: &mut Document, source: &str) -> Result<(Film, f64), String> {
    let context = FilmContext::from_document(doc);
    let scale = context.unit_scale;
    FilmHost::new(HostOptions::default())
        .run(
            source,
            "guide.rhai",
            &film_name("films/guide.rhai"),
            context,
        )
        .map(|run| (run.film, scale))
        .map_err(|e| e.render())
}

#[test]
fn the_guide_written_with_insert_and_screw_in_moves_every_screw_as_its_explicit_tracks() {
    let mut doc = angled_fitting();
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let explicit =
        std::fs::read_to_string(dir.join("tests/fixtures/angled_fitting.guide.explicit.film.rhai"))
            .expect("explicit guide");
    let insert = std::fs::read_to_string(dir.join("examples/angled_fitting.guide.film.rhai"))
        .expect("guide");
    let (a, scale) = film(&mut doc, &explicit).expect("explicit film");
    let (b, _) = film(&mut doc, &insert).expect("insert film");
    assert!((a.duration - b.duration).abs() < 1e-12);
    let frames = (a.duration * f64::from(a.fps)).ceil() as u32;
    let mut worst: f64 = 0.0;
    for screw in ["angled_screw", "mount_screw_a", "mount_screw_b"] {
        for frame in 0..=frames {
            let t = f64::from(frame) / f64::from(a.fps);
            let (pa, pb) = (
                occurrence_motion(&a, screw, t, scale).expect("motion"),
                occurrence_motion(&b, screw, t, scale).expect("motion"),
            );
            let (ma, mb) = match (pa, pb) {
                (Some(x), Some(y)) => (x.transform().matrix, y.transform().matrix),
                (None, None) => continue,
                other => panic!("{screw} at {t}: moved in one film only: {other:?}"),
            };
            for r in 0..3 {
                for c in 0..4 {
                    worst = worst.max((ma[r][c] - mb[r][c]).abs());
                }
            }
        }
    }
    assert!(worst < 1e-9, "poses differ by {worst:e}");
}

#[test]
fn an_occurrence_with_no_axis_declines_naming_it() {
    let mut doc = angled_fitting();
    let error = film(
        &mut doc,
        r#"let film = storyboard();
film.step("in", "", 1.0, |s| { s.insert("bracket_1", #{ distance: 5.0 }); });
film"#,
    )
    .expect_err("the bracket is no fastener");
    assert!(
        error.contains("bracket_1") && error.contains("not a fastener"),
        "{error}"
    );
}
