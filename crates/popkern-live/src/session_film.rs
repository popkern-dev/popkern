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
    /// against `shown` — the context of the model the page shows (a
    /// pre-built envelope's, or a live one's) — or, when the page sends
    /// none, the session's built document. The page's context wins: a
    /// pre-built model is not in the session, which still holds the
    /// last model built live, and a film read that one's targets (a
    /// hinged box's film saw the hinged case's `box` and `lid`; PK-869).
    pub fn film(
        &mut self,
        source: &str,
        name: &str,
        shown: Option<FilmContext>,
    ) -> Result<FilmReply, Outcome> {
        let options = self.options();
        let context = match (shown, self.document_mut()) {
            (Some(context), _) => context,
            (None, Some(document)) => FilmContext::from_document(document),
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
