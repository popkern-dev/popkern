//! Presentation films (docs/architecture/presentation-films.md): run a
//! film script against the session's built document. The document is
//! read, never rebuilt or changed, so a film reruns in milliseconds.

use popkern_rhai::film::{Film, FilmContext, FilmHost, film_name};
use popkern_rhai::outcome::from_run_error;
use serde::Serialize;

use crate::session::Session;
use popkern_rhai::Outcome;

/// A film run's reply.
#[derive(Debug, Clone, Serialize)]
pub struct FilmReply {
    /// The storyboard.
    pub film: Film,
    /// Lines the film script printed.
    pub printed: Vec<String>,
    /// What the film could name (occurrences, bodies, joints, unit).
    pub context: FilmContext,
}

impl Session {
    /// Run the film script `source` (file `name`: `films/guide.rhai`)
    /// against the built document, or against `fallback` (the context
    /// a pre-built envelope carried) when the session holds none.
    pub fn film(
        &mut self,
        source: &str,
        name: &str,
        fallback: Option<FilmContext>,
    ) -> Result<FilmReply, Outcome> {
        let options = self.options();
        let context = match (self.document_mut(), fallback) {
            (Some(document), _) => FilmContext::from_document(document),
            (None, Some(context)) => context,
            (None, None) => {
                return Err(Outcome::GateFailure {
                    step: "film".into(),
                    detail: "build the model first: a film plays the built document".into(),
                });
            }
        };
        let host = FilmHost::new(popkern_rhai::HostOptions {
            track_sources: false,
            ..options
        });
        let script = name.rsplit('/').next().unwrap_or(name);
        let run = host
            .run(source, script, &film_name(name), context.clone())
            .map_err(|error| from_run_error(&error))?;
        Ok(FilmReply {
            film: run.film,
            printed: run.printed,
            context,
        })
    }
}
