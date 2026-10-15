//! Presentation films (docs/architecture/presentation-films.md): a film
//! script reads a snapshot of the built document and writes a storyboard
//! the page plays and exports. Nothing here touches the document.

mod actions;
mod api;
mod context;
pub mod model;
mod register;

use std::cell::RefCell;
use std::rc::Rc;

use rhai::{Dynamic, Scope};

pub use api::{ScriptFilm, ScriptFilmDocument, ScriptStep};
pub use context::FilmContext;
pub use model::{CameraKey, CameraTo, Channel, Ease, Film, FilmStep, SectionKey, Track};

use crate::engine::build_engine;
use crate::host::{HostOptions, ScriptRunError};

/// What a film run produced.
#[derive(Debug, Clone)]
pub struct FilmRun {
    /// The storyboard.
    pub film: Film,
    /// Lines the film script printed.
    pub printed: Vec<String>,
}

/// Runs film scripts against a document snapshot.
pub struct FilmHost {
    options: HostOptions,
}

impl FilmHost {
    /// A host with the given budgets (operations and wall clock).
    pub fn new(options: HostOptions) -> FilmHost {
        FilmHost { options }
    }

    /// Run the film script `source` (named `script` in diagnostics, its
    /// film called `name`) against `context`.
    pub fn run(
        &self,
        source: &str,
        script: &str,
        name: &str,
        context: FilmContext,
    ) -> Result<FilmRun, ScriptRunError> {
        let printed: Rc<RefCell<Vec<String>>> = Rc::new(RefCell::new(Vec::new()));
        let sink = Rc::clone(&printed);
        let mut engine = build_engine(move |line| sink.borrow_mut().push(line.to_string()));
        if let Some(max) = self.options.max_operations {
            engine.set_max_operations(max);
        }
        register::register(&mut engine, Rc::new(context), name);
        let watch = popkern_core::time::Stopwatch::start();
        let timeout = self.options.timeout;
        engine.on_progress(move |_| match timeout {
            Some(timeout) if watch.elapsed().is_some_and(|elapsed| elapsed > timeout) => Some(
                Dynamic::from(format!("wall-clock budget of {timeout:?} exceeded")),
            ),
            _ => None,
        });
        let ast = engine.compile(source).map_err(|e| {
            let error: rhai::EvalAltResult = e.into();
            ScriptRunError::from_eval(&error, script)
        })?;
        let mut scope = Scope::new();
        let returned = engine
            .eval_ast_with_scope::<Dynamic>(&mut scope, &ast)
            .map_err(|e| ScriptRunError::from_eval(&e, script))?;
        // The returned storyboard, else the last one left in scope.
        let film = returned
            .try_cast::<ScriptFilm>()
            .or_else(|| {
                scope
                    .iter()
                    .filter_map(|(_, _, v)| v.try_cast::<ScriptFilm>())
                    .last()
            })
            .map(|f| f.film());
        let Some(film) = film else {
            return Err(ScriptRunError::argument(
                script,
                "a film script returns its storyboard: end it with `film` (let film = storyboard(...))",
            ));
        };
        if film.steps.is_empty() {
            return Err(ScriptRunError::argument(
                script,
                "the storyboard has no steps: add film.step(title, seconds, |s| { ... })",
            ));
        }
        let printed = printed.borrow().clone();
        Ok(FilmRun { film, printed })
    }
}

/// A film's name from its file name: `films/assembly_guide.rhai` →
/// `assembly_guide`; `angled_fitting.guide.film.rhai` → `guide`;
/// `angled_fitting.film.rhai` → `angled_fitting`.
pub fn film_name(path: &str) -> String {
    let file = path.rsplit('/').next().unwrap_or(path);
    let stem = file.strip_suffix(".rhai").unwrap_or(file);
    match stem.strip_suffix(".film") {
        Some(rest) => rest.rsplit('.').next().unwrap_or(rest).to_string(),
        None => stem.to_string(),
    }
}
