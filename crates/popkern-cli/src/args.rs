//! Command-line parsing (hand-rolled: the surface is small and the
//! workspace takes no new dependency families for it).

use std::path::PathBuf;
use std::time::Duration;

/// Usage text.
pub const USAGE: &str = "\
usage: popkern run <script.rhai> [options]

  --step <path>        export the active body as STEP
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

exit codes: 0 served, 2 script error, 3 kernel decline, 4 gate or export failure";

/// The subcommand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Command {
    /// `popkern run`.
    Run,
    /// `popkern help` / `--help`.
    Help,
}

/// Parsed arguments.
#[derive(Debug, Clone)]
pub struct Cli {
    /// The subcommand.
    pub command: Command,
    /// The script path.
    pub script: PathBuf,
    /// STEP export path.
    pub step: Option<PathBuf>,
    /// PNG render path.
    pub png: Option<PathBuf>,
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
}

impl Default for Cli {
    fn default() -> Self {
        Cli {
            command: Command::Run,
            script: PathBuf::new(),
            step: None,
            png: None,
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

/// Parse the arguments after the program name.
pub fn parse_args(args: &[String]) -> Result<Cli, UsageError> {
    let mut cli = Cli::default();
    let mut iter = args.iter();
    match iter.next().map(String::as_str) {
        Some("run") => {}
        Some("help") | Some("--help") | Some("-h") => {
            cli.command = Command::Help;
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
            "--pkn" => cli.pkn = Some(PathBuf::from(value(arg, &mut iter)?)),
            "--json" => cli.json = true,
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
            other if other.starts_with("--") => {
                return Err(UsageError(format!("unknown option '{other}'")));
            }
            path => {
                if script.replace(PathBuf::from(path)).is_some() {
                    return Err(UsageError("only one script may be given".into()));
                }
            }
        }
    }
    cli.script = script.ok_or_else(|| UsageError("missing script path".into()))?;
    Ok(cli)
}
