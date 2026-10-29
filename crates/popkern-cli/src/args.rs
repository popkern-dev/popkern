//! Command-line parsing (hand-rolled: the surface is small and the
//! workspace takes no new dependency families for it).

use std::path::PathBuf;
use std::time::Duration;

/// Usage text.
pub const USAGE: &str = "\
usage: popkern run <script.rhai> [options]
       popkern build <project-dir> [--pkn <path>] [--step <path>] [--stl <path>] [--json]
                     [--strict] [--rebuild] [--max-ops <n>] [--timeout <s>]
       popkern new <project-dir> [--name <text>]
       popkern serve <projects-dir> [--page <dir>] [--port <n>] [--bind <addr>]
                     [--compute native|browser]
       popkern film <project-dir | script.rhai> [--film <name>] [--json]
                     [--out <file> [--preset <name>] [--at <s>] [--chrome <path>] [--page <dir>]]
       popkern --version

popkern build runs a project's migrations (migrations/NNN_name.rhai) in
order into one document. With --pkn naming an existing document, only
the migrations it has not applied run; a document whose applied history
differs from the folder is rebuilt (--strict: refused, exit 6).
--rebuild always rebuilds from the folder.

popkern serve runs popkern live self-hosted: projects are the folders
under <projects-dir>, and the kernel runs natively (--compute browser
keeps it in the page). --page names the page folder (default: the
installed page, share/popkern/live beside the binary's bin/ folder, else
web/live built with tools/build-live.sh). Local only: --bind must be a loopback
address (default 127.0.0.1, port 7690). The same server answers MCP at
/mcp (JSON-RPC over HTTP): agents list, edit, build and export the same
projects, e.g.  claude mcp add --transport http popkern http://127.0.0.1:7690/mcp

run options:

  --step <path>        export the model as STEP (every body)
  --stl <path>         export every visible body as one binary STL (also build)
  --stl-chord <len>    STL chord height, document units (default 0.005)
  --stl-angle <deg>    STL facet angle bound (default 12)
  --png <path>         render the active body to PNG
  --pkn <path>         save the document
  --json               print the report as JSON instead of text
  --param <key=value>  set a declared parameter after the run and replay
  --tolerance <t>      linear tolerance for the document (overrides the script)
  --max-ops <n>        stop the script after n interpreter operations
  --timeout <seconds>  stop the script after a wall-clock budget
  --watch              rebuild on save (Ctrl-C to stop)
  --width <px>         render width (default 800)
  --height <px>        render height (default 600)
  --no-edges           render without edge lines
  --zebra[=mesh|surface]
                       render --png as reflection-line zebra stripes, from
                       the display mesh's normals (mesh, the default) or the
                       surfaces' own (surface), and report how far the two
                       disagree
  --stripes <n>        zebra stripe count (default 16)
  --facet-angle <deg>  zebra mesh facet bound (default 12, the display bound)

popkern film runs a project's films (films/<name>.rhai) or a script's
(stem.<name>.film.rhai beside it) against the built document: --json
prints the storyboards. --out renders one film to a file with headless
Chrome (--chrome or POPKERN_CHROME): presets mp4-1080p (default),
mp4-vertical, mp4-720p, webm, gif, png (one frame, at --at seconds).

exit codes: 0 served, 2 script error, 3 kernel decline, 4 gate or export failure,
            5 project structure error, 6 diverged history under --strict";

/// The subcommand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Command {
    /// `popkern run`.
    Run,
    /// `popkern build <dir>`.
    Build,
    /// `popkern new <dir>`.
    New,
    /// `popkern serve <dir>`.
    Serve,
    /// `popkern film <dir | script>`.
    Film,
    /// `popkern help` / `--help`.
    Help,
    /// `popkern --version` / `version`.
    Version,
}

/// Parsed arguments.
#[derive(Debug, Clone)]
pub struct Cli {
    /// The subcommand.
    pub command: Command,
    /// The script path (`run`) or project folder (`build`, `new`).
    pub script: PathBuf,
    /// STEP export path.
    pub step: Option<PathBuf>,
    /// PNG render path.
    pub png: Option<PathBuf>,
    /// STL export path.
    pub stl: Option<PathBuf>,
    /// STL chord height override (document units).
    pub stl_chord: Option<f64>,
    /// STL facet angle bound override (degrees).
    pub stl_angle: Option<f64>,
    /// Document save path.
    pub pkn: Option<PathBuf>,
    /// JSON report.
    pub json: bool,
    /// Parameter overrides applied after the run.
    pub params: Vec<(String, f64)>,
    /// Linear tolerance override.
    pub tolerance: Option<f64>,
    /// Operation budget.
    pub max_ops: Option<u64>,
    /// Wall-clock budget.
    pub timeout: Option<Duration>,
    /// Rebuild on save.
    pub watch: bool,
    /// Render size.
    pub width: u32,
    /// Render size.
    pub height: u32,
    /// Draw edges in the render.
    pub edges: bool,
    /// Render the PNG as a zebra view in this normal source.
    pub zebra: Option<popkern_viz::zebra::ZebraNormals>,
    /// Zebra stripe count.
    pub stripes: u32,
    /// Zebra mesh angular facet bound (degrees).
    pub facet_angle_deg: f64,
    /// `build`: refuse a saved document whose history diverged.
    pub strict: bool,
    /// `build`: rebuild from the folder whatever the saved document says.
    pub rebuild: bool,
    /// `new`: the project's display name.
    pub name: Option<String>,
    /// `serve`: the page folder.
    pub page: Option<PathBuf>,
    /// `serve`: the port.
    pub port: u16,
    /// `serve`: the address to bind.
    pub bind: String,
    /// `serve`: where the kernel runs (`native` or `browser`).
    pub compute: String,
    /// `popkern film --film`: the film to run (default: every film, or
    /// the only one for `--out`).
    pub film: Option<String>,
    /// `popkern film --out`: the file to render the film to.
    pub out: Option<PathBuf>,
    /// `popkern film --preset`: the export preset.
    pub preset: String,
    /// `popkern film --chrome`: the headless Chrome binary.
    pub chrome: Option<PathBuf>,
    /// `popkern film --at`: the time of a `png` frame, seconds.
    pub at: f64,
}

impl Default for Cli {
    fn default() -> Self {
        Cli {
            command: Command::Run,
            script: PathBuf::new(),
            step: None,
            png: None,
            stl: None,
            stl_chord: None,
            stl_angle: None,
            pkn: None,
            json: false,
            params: Vec::new(),
            tolerance: None,
            max_ops: None,
            timeout: None,
            watch: false,
            width: 800,
            height: 600,
            edges: true,
            zebra: None,
            stripes: 16,
            facet_angle_deg: popkern_mesh::DISPLAY_MAX_FACET_ANGLE.to_degrees(),
            strict: false,
            rebuild: false,
            name: None,
            page: None,
            port: 7690,
            bind: "127.0.0.1".into(),
            compute: "native".into(),
            film: None,
            out: None,
            preset: "mp4-1080p".into(),
            chrome: None,
            at: 0.0,
        }
    }
}

/// A usage error.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("{0}")]
pub struct UsageError(pub String);

fn value<'a>(
    flag: &str,
    iter: &mut impl Iterator<Item = &'a String>,
) -> Result<&'a String, UsageError> {
    iter.next()
        .ok_or_else(|| UsageError(format!("{flag} needs a value")))
}

fn number<T: std::str::FromStr>(flag: &str, text: &str) -> Result<T, UsageError> {
    text.parse()
        .map_err(|_| UsageError(format!("{flag}: '{text}' is not a number")))
}

fn positive(flag: &str, text: &str) -> Result<f64, UsageError> {
    let v: f64 = number(flag, text)?;
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(UsageError(format!("{flag}: '{text}' must be positive")))
    }
}

/// Parse the arguments after the program name.
pub fn parse_args(args: &[String]) -> Result<Cli, UsageError> {
    let mut cli = Cli::default();
    let mut iter = args.iter();
    match iter.next().map(String::as_str) {
        Some("run") => {}
        Some("build") => cli.command = Command::Build,
        Some("new") => cli.command = Command::New,
        Some("serve") => cli.command = Command::Serve,
        Some("film") => cli.command = Command::Film,
        Some("help") | Some("--help") | Some("-h") => {
            cli.command = Command::Help;
            return Ok(cli);
        }
        Some("version") | Some("--version") | Some("-V") => {
            cli.command = Command::Version;
            return Ok(cli);
        }
        Some(other) => return Err(UsageError(format!("unknown command '{other}'"))),
        None => return Err(UsageError("missing command".into())),
    }
    let mut script = None;
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--step" => cli.step = Some(PathBuf::from(value(arg, &mut iter)?)),
            "--png" => cli.png = Some(PathBuf::from(value(arg, &mut iter)?)),
            "--stl" => cli.stl = Some(PathBuf::from(value(arg, &mut iter)?)),
            "--stl-chord" => cli.stl_chord = Some(positive(arg, value(arg, &mut iter)?)?),
            "--stl-angle" => cli.stl_angle = Some(positive(arg, value(arg, &mut iter)?)?),
            "--pkn" => cli.pkn = Some(PathBuf::from(value(arg, &mut iter)?)),
            "--json" => cli.json = true,
            "--strict" => cli.strict = true,
            "--rebuild" => cli.rebuild = true,
            "--name" => cli.name = Some(value(arg, &mut iter)?.clone()),
            "--page" => cli.page = Some(PathBuf::from(value(arg, &mut iter)?)),
            "--port" => cli.port = number(arg, value(arg, &mut iter)?)?,
            "--bind" => cli.bind = value(arg, &mut iter)?.clone(),
            "--compute" => {
                let v = value(arg, &mut iter)?;
                if v != "native" && v != "browser" {
                    return Err(UsageError(format!(
                        "--compute is native or browser, not '{v}'"
                    )));
                }
                cli.compute = v.clone();
            }
            "--film" => cli.film = Some(value(arg, &mut iter)?.clone()),
            "--out" => cli.out = Some(PathBuf::from(value(arg, &mut iter)?)),
            "--preset" => cli.preset = value(arg, &mut iter)?.clone(),
            "--chrome" => cli.chrome = Some(PathBuf::from(value(arg, &mut iter)?)),
            "--at" => cli.at = number(arg, value(arg, &mut iter)?)?,
            "--watch" => cli.watch = true,
            "--no-edges" => cli.edges = false,
            "--param" => {
                let text = value(arg, &mut iter)?;
                let (key, v) = text.split_once('=').ok_or_else(|| {
                    UsageError(format!("--param expects key=value, got '{text}'"))
                })?;
                cli.params
                    .push((key.trim().to_string(), number("--param", v.trim())?));
            }
            "--tolerance" => cli.tolerance = Some(number(arg, value(arg, &mut iter)?)?),
            "--max-ops" => cli.max_ops = Some(number(arg, value(arg, &mut iter)?)?),
            "--timeout" => {
                let seconds: f64 = number(arg, value(arg, &mut iter)?)?;
                cli.timeout = Some(Duration::from_secs_f64(seconds.max(0.0)));
            }
            "--width" => cli.width = number(arg, value(arg, &mut iter)?)?,
            "--height" => cli.height = number(arg, value(arg, &mut iter)?)?,
            "--zebra" | "--zebra=mesh" => cli.zebra = Some(popkern_viz::zebra::ZebraNormals::Mesh),
            "--zebra=surface" => cli.zebra = Some(popkern_viz::zebra::ZebraNormals::Surface),
            "--stripes" => cli.stripes = number(arg, value(arg, &mut iter)?)?,
            "--facet-angle" => cli.facet_angle_deg = number(arg, value(arg, &mut iter)?)?,
            other if other.starts_with("--") => {
                return Err(UsageError(format!("unknown option '{other}'")));
            }
            path => {
                if script.replace(PathBuf::from(path)).is_some() {
                    return Err(UsageError("only one script or folder may be given".into()));
                }
            }
        }
    }
    let what = match cli.command {
        Command::Build | Command::New => "missing project folder",
        Command::Serve => "missing projects folder",
        Command::Film => "missing project folder or script",
        _ => "missing script path",
    };
    cli.script = script.ok_or_else(|| UsageError(what.into()))?;
    Ok(cli)
}
