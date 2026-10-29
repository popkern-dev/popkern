//! Native film output and still renders (cli-render-export.md §2, §5):
//! `popkern film --out` without a browser, and `popkern render`.

use popkern_rhai::film::state::frame_times;

use crate::args::Cli;
use crate::film::{FilmRun, load, run_films};
use crate::film_caption::Fonts;
use crate::film_encode::{Format, GifWriter, Mp4Writer, png, preset};
use crate::film_native::{NativeFilm, still_script};

/// The caption fonts a film needs: none when it has no caption; a
/// warning (and captionless frames) when no font is found.
fn fonts_for(cli: &Cli, film: &popkern_rhai::film::Film) -> Option<Fonts> {
    if film
        .steps
        .iter()
        .all(|s| s.title.is_empty() || s.kind == "hold")
    {
        return None;
    }
    match Fonts::find(cli.font.as_deref()) {
        Ok(fonts) => {
            eprintln!("captions in {}", fonts.source.display());
            Some(fonts)
        }
        Err(message) => {
            eprintln!("popkern: {message}; captions are left out");
            None
        }
    }
}

/// Render film `index` of `run` with `cli.preset`: the file's bytes.
pub fn render_native(cli: &Cli, run: &mut FilmRun, index: usize) -> Result<Vec<u8>, String> {
    let preset = preset(&cli.preset).ok_or_else(|| format!("no preset {}", cli.preset))?;
    let film = run
        .films
        .get(index)
        .map(|(_, f)| f.clone())
        .ok_or("no such film")?;
    let [mut width, mut height] = preset.frame_size(film.size);
    if preset.format == Format::Mp4 {
        width &= !1;
        height &= !1;
    }
    let fps = preset.fps.unwrap_or(film.fps).max(1);
    let fonts = fonts_for(cli, &film);
    let envelope = run.envelope.clone();
    let mut native = NativeFilm::new(&mut run.session, envelope, film, fonts)?;
    if preset.format == Format::Png {
        let frame = native.frame(cli.at, width, height)?;
        return png(&frame, width, height);
    }
    let times = frame_times(native.film(), fps as f64);
    let started = popkern_core::time::Stopwatch::start();
    eprintln!(
        "rendering {} frame(s) at {width}x{height}, {fps} fps",
        times.len()
    );
    let report = |k: usize| {
        if (k + 1).is_multiple_of(fps as usize) || k + 1 == times.len() {
            eprint!(
                "\r  {}/{} frames ({:.0} s)",
                k + 1,
                times.len(),
                started.elapsed_ms().unwrap_or(0.0) / 1000.0
            );
        }
    };
    let bytes = match preset.format {
        Format::Gif => {
            let mut gif = GifWriter::new(width, height, fps)?;
            for (k, &t) in times.iter().enumerate() {
                gif.frame(&native.frame(t, width, height)?)?;
                report(k);
            }
            gif.finish()?
        }
        Format::Mp4 => {
            let mut mp4 = Mp4Writer::new(width, height, fps, preset.bitrate)?;
            for (k, &t) in times.iter().enumerate() {
                mp4.frame(&native.frame(t, width, height)?)?;
                report(k);
            }
            mp4.finish()?
        }
        Format::Webm => {
            return Err(
                "WebM is rendered by the browser: add --browser, or use an mp4 \
                        or gif preset"
                    .into(),
            );
        }
        Format::Png => Vec::new(),
    };
    eprintln!();
    Ok(bytes)
}

/// `popkern render`: one image of a script or project. Exit code.
pub fn render_command(cli: &Cli) -> u8 {
    let Some(out) = &cli.out else {
        eprintln!("popkern render: --out <file.png> is required");
        return 2;
    };
    let mut source = match load(&cli.script) {
        Ok(source) => source,
        Err(message) => {
            eprintln!("popkern render: {message}");
            return 5;
        }
    };
    let mut cli = cli.clone();
    if cli.film.is_none() {
        // A still: a one-frame film holding the requested view.
        let path = "films/still.rhai".to_string();
        let text = still_script(
            cli.look.as_deref().unwrap_or("studio"),
            cli.background.as_deref().unwrap_or("#ffffff"),
            cli.azimuth,
            cli.elevation,
            cli.zoom,
        );
        source.files.push((path.clone(), text.into_bytes()));
        source.films = vec![("still".into(), path)];
        cli.film = Some("still".into());
        cli.at = 0.0;
    }
    let mut run = match run_films(&cli, &source) {
        Ok(run) => run,
        Err((code, message)) => {
            eprintln!("popkern render: {message}");
            return code;
        }
    };
    let Some(film) = run.films.first().map(|(_, f)| f.clone()) else {
        eprintln!("popkern render: no film");
        return 2;
    };
    let fonts = fonts_for(&cli, &film);
    let envelope = run.envelope.clone();
    let rendered = NativeFilm::new(&mut run.session, envelope, film, fonts)
        .and_then(|mut native| native.frame(cli.at, cli.width, cli.height))
        .and_then(|frame| png(&frame, cli.width, cli.height))
        .and_then(|bytes| {
            crate::render::create_parent(out)?;
            std::fs::write(out, &bytes).map_err(|e| format!("{}: {e}", out.display()))
        });
    match rendered {
        Ok(()) => {
            println!("wrote {} ({}x{})", out.display(), cli.width, cli.height);
            0
        }
        Err(message) => {
            eprintln!("popkern render: {message}");
            4
        }
    }
}
