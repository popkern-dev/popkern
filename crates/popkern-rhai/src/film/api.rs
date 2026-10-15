//! The film scripting surface (docs/architecture/presentation-films.md
//! §4): `document()`, `storyboard(options)`, the step builders and every
//! action, each validated against the built document at its own line.

use std::cell::RefCell;
use std::rc::Rc;

use rhai::{Dynamic, EvalAltResult, FnPtr, NativeCallContext};

use super::context::FilmContext;
use super::model::{CameraKey, Channel, Ease, Film, FilmStep, Track};
use crate::args::{self, Options};
use crate::error::argument_error;

/// Longest step, seconds.
const MAX_STEP_SECONDS: f64 = 600.0;
/// Largest frame side, pixels.
const MAX_SIDE: u32 = 3840;

/// The film being written, shared by the storyboard value and its steps.
#[derive(Debug)]
pub(crate) struct FilmState {
    pub(crate) film: Film,
    pub(crate) context: Rc<FilmContext>,
}

/// `storyboard(...)`'s value.
#[derive(Debug, Clone)]
pub struct ScriptFilm(pub(crate) Rc<RefCell<FilmState>>);

impl ScriptFilm {
    /// The film as written so far.
    pub fn film(&self) -> Film {
        self.0.borrow().film.clone()
    }
}

/// The `s` a step's closure receives.
#[derive(Debug, Clone)]
pub struct ScriptStep {
    state: Rc<RefCell<FilmState>>,
    pub(super) start: f64,
    pub(super) duration: f64,
}

/// `document()`'s value: the read-only snapshot.
#[derive(Debug, Clone)]
pub struct ScriptFilmDocument(pub(super) Rc<FilmContext>);

fn line(ctx: &NativeCallContext) -> Option<usize> {
    ctx.call_position().line()
}

pub(super) fn seconds(argument: &str, value: &Dynamic) -> Result<f64, Box<EvalAltResult>> {
    let s = args::float(argument, value)?;
    if !(s.is_finite() && (0.0..=MAX_STEP_SECONDS).contains(&s)) {
        return Err(argument_error(
            argument,
            format!("seconds between 0 and {MAX_STEP_SECONDS}, got {s}"),
        ));
    }
    Ok(s)
}

pub(super) fn ease(argument: &str, value: &Dynamic) -> Result<Ease, Box<EvalAltResult>> {
    match args::text(argument, value)?.as_str() {
        "linear" => Ok(Ease::Linear),
        "in" => Ok(Ease::In),
        "out" => Ok(Ease::Out),
        "inout" | "in_out" => Ok(Ease::InOut),
        other => Err(argument_error(
            argument,
            format!("\"linear\", \"in\", \"out\" or \"inout\", got \"{other}\""),
        )),
    }
}

pub(super) fn target(
    context: &FilmContext,
    argument: &str,
    value: &Dynamic,
) -> Result<String, Box<EvalAltResult>> {
    let name = args::text(argument, value)?;
    if context.selects(&name) {
        return Ok(name);
    }
    Err(argument_error(
        argument,
        format!(
            "the document has no '{name}'; it has {}",
            context.known_targets()
        ),
    ))
}

pub(super) fn length3(
    context: &FilmContext,
    argument: &str,
    value: &Dynamic,
) -> Result<[f64; 3], Box<EvalAltResult>> {
    let v = args::vec3(argument, value)?;
    Ok(v.map(|c| c * context.unit_scale))
}

pub(super) fn unit3(argument: &str, value: &Dynamic) -> Result<[f64; 3], Box<EvalAltResult>> {
    let v = args::vec3(argument, value)?;
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if !(n.is_finite() && n > 1e-12) {
        return Err(argument_error(argument, "a non-zero direction"));
    }
    Ok(v.map(|c| c / n))
}

const STORYBOARD_OPTIONS: &[&str] = &["size", "fps", "look", "background", "camera"];
const CAMERA_OPTIONS: &[&str] = &["azimuth", "elevation", "zoom", "target"];

fn camera_key(
    context: &FilmContext,
    base: &CameraKey,
    options: &Options,
) -> Result<CameraKey, Box<EvalAltResult>> {
    let zoom = options.with("zoom", args::float)?.unwrap_or(base.zoom);
    if !(zoom.is_finite() && zoom > 0.0) {
        return Err(argument_error(
            "zoom",
            format!("a positive factor, got {zoom}"),
        ));
    }
    let elevation = options
        .with("elevation", args::float)?
        .unwrap_or(base.elevation);
    if !(-89.0..=89.0).contains(&elevation) {
        return Err(argument_error(
            "elevation",
            format!("degrees in -89..89, got {elevation}"),
        ));
    }
    Ok(CameraKey {
        azimuth: options
            .with("azimuth", args::float)?
            .unwrap_or(base.azimuth),
        elevation,
        zoom,
        target: match options.get("target") {
            Some(v) => target(context, "target", v)?,
            None => base.target.clone(),
        },
    })
}

fn colour(argument: &str, value: &Dynamic) -> Result<String, Box<EvalAltResult>> {
    let text = args::text(argument, value)?;
    let ok = text.len() == 7
        && text.starts_with('#')
        && text[1..].chars().all(|c| c.is_ascii_hexdigit());
    if !ok {
        return Err(argument_error(
            argument,
            format!("a colour \"#rrggbb\", got \"{text}\""),
        ));
    }
    Ok(text.to_ascii_lowercase())
}

pub(super) fn storyboard(
    context: Rc<FilmContext>,
    name: &str,
    options: &Dynamic,
) -> Result<ScriptFilm, Box<EvalAltResult>> {
    let options = Options::parse("storyboard", options, STORYBOARD_OPTIONS)?;
    let size = match options.get("size") {
        Some(v) => {
            let list = v
                .clone()
                .into_typed_array::<Dynamic>()
                .map_err(|_| argument_error("size", "[width, height] in pixels"))?;
            if list.len() != 2 {
                return Err(argument_error("size", "[width, height] in pixels"));
            }
            let side = |k: usize| -> Result<u32, Box<EvalAltResult>> {
                let s = args::float("size", &list[k])?;
                if !(s.is_finite() && s >= 16.0 && s <= MAX_SIDE as f64) {
                    return Err(argument_error(
                        "size",
                        format!("sides between 16 and {MAX_SIDE} pixels, got {s}"),
                    ));
                }
                Ok((s.round() as u32) & !1)
            };
            [side(0)?, side(1)?]
        }
        None => [1920, 1080],
    };
    let fps = match options.with("fps", args::float)? {
        Some(f) if (1.0..=120.0).contains(&f) => f.round() as u32,
        Some(f) => {
            return Err(argument_error(
                "fps",
                format!("1..120 frames per second, got {f}"),
            ));
        }
        None => 30,
    };
    let look = match options.with("look", args::text)? {
        Some(l) if ["toon", "studio", "cad"].contains(&l.as_str()) => l,
        Some(l) => {
            return Err(argument_error(
                "look",
                format!("\"toon\", \"studio\" or \"cad\", got \"{l}\""),
            ));
        }
        None => "toon".into(),
    };
    let background = options
        .with("background", colour)?
        .unwrap_or_else(|| "#ffffff".into());
    let base = CameraKey {
        azimuth: -35.0,
        elevation: 25.0,
        zoom: 1.0,
        target: "all".into(),
    };
    let camera = match options.get("camera") {
        Some(v) => {
            let o = Options::parse("camera", v, CAMERA_OPTIONS)?;
            camera_key(&context, &base, &o)?
        }
        None => base,
    };
    let film = Film {
        name: name.to_string(),
        size,
        fps,
        look,
        background,
        duration: 0.0,
        camera,
        steps: Vec::new(),
        tracks: Vec::new(),
    };
    Ok(ScriptFilm(Rc::new(RefCell::new(FilmState {
        film,
        context,
    }))))
}

/// Append a step and run its closure.
pub(super) fn add_step(
    ctx: &NativeCallContext,
    film: &mut ScriptFilm,
    kind: &'static str,
    title: &str,
    subtitle: &str,
    duration: &Dynamic,
    body: Option<FnPtr>,
) -> Result<ScriptFilm, Box<EvalAltResult>> {
    let duration = seconds("seconds", duration)?;
    if duration <= 0.0 {
        return Err(argument_error(
            "seconds",
            "a step lasts longer than 0 seconds",
        ));
    }
    let start = {
        let mut state = film.0.borrow_mut();
        let start = state.film.duration;
        state.film.steps.push(FilmStep {
            title: title.to_string(),
            subtitle: subtitle.to_string(),
            start,
            duration,
            kind,
            line: line(ctx),
        });
        state.film.duration = start + duration;
        start
    };
    if let Some(body) = body {
        let step = ScriptStep {
            state: Rc::clone(&film.0),
            start,
            duration,
        };
        let _: Dynamic = body.call_within_context(ctx, (step,))?;
    }
    Ok(film.clone())
}

impl ScriptStep {
    pub(super) fn context(&self) -> Rc<FilmContext> {
        Rc::clone(&self.state.borrow().context)
    }

    /// Timing from the action's options: `at` into the step, `duration`
    /// (default to the step's end), `ease`.
    pub(super) fn timing(&self, options: &Options) -> Result<(f64, f64, Ease), Box<EvalAltResult>> {
        let at = options.with("at", seconds)?.unwrap_or(0.0);
        let rest = (self.duration - at).max(0.0);
        let duration = options.with("duration", seconds)?.unwrap_or(rest);
        let ease = options.with("ease", ease)?.unwrap_or(Ease::InOut);
        Ok((self.start + at, duration, ease))
    }

    pub(super) fn push(&self, ctx: &NativeCallContext, timing: (f64, f64, Ease), channel: Channel) {
        let (start, duration, ease) = timing;
        self.state.borrow_mut().film.tracks.push(Track {
            start,
            duration,
            ease,
            line: line(ctx),
            channel,
        });
    }
}
