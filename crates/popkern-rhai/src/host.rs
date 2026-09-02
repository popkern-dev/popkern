//! The script host: builds the locked-down engine, runs a script,
//! and hands back the model the script authored plus everything it
//! printed. Budgets (operation count, wall clock) are host options so
//! a CLI, a daemon, and a browser page bound scripts the same way.

use std::cell::RefCell;
use std::rc::Rc;
use std::time::Duration;

use popkern_core::time::Stopwatch;
use rhai::{AST, Dynamic, Engine, EvalAltResult, Position, Scope};
use serde_json::Value;

use crate::engine::build_engine;
use crate::model_api::ScriptModel;

/// Limits a host places on a script run.
#[derive(Debug, Clone, Copy, Default)]
pub struct HostOptions {
    /// Maximum Rhai operations before the run is terminated
    /// (`None` = unlimited).
    pub max_operations: Option<u64>,
    /// Wall-clock budget; checked at every Rhai progress tick and
    /// honoured only where a monotonic clock exists (`Stopwatch`).
    pub timeout: Option<Duration>,
}

/// What a finished run produced.
pub struct RunOutcome {
    /// The model the script authored: the script's return value when
    /// it is a model, else the last `Model` variable in scope.
    pub model: Option<popkern_script::Model>,
    /// Lines the script printed, in order.
    pub printed: Vec<String>,
    /// Wall time of the run when a clock was available.
    pub elapsed: Option<Duration>,
}

/// A failed run: the §10 payload plus where in the script it was
/// raised. `kind` is the payload class (`BlendDeclined`,
/// `ArgumentError`, …), `SyntaxError` for a parse failure, and
/// `BudgetExceeded` when a host limit terminated the run.
#[derive(Debug, Clone, thiserror::Error)]
pub struct ScriptRunError {
    /// The §10 class.
    pub kind: String,
    /// The `Display` text of the underlying error.
    pub message: String,
    /// Script name the host passed to `run`.
    pub script: String,
    /// 1-based line of the raising call, when known.
    pub line: Option<usize>,
    /// 1-based column of the raising call, when known.
    pub column: Option<usize>,
    /// The full payload map (the `catch` value) as JSON.
    pub payload: Box<Value>,
}

impl std::fmt::Debug for RunOutcome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RunOutcome")
            .field("model", &self.model.is_some())
            .field("printed", &self.printed)
            .field("elapsed", &self.elapsed)
            .finish()
    }
}

impl std::fmt::Display for ScriptRunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.kind)?;
        if let Some(line) = self.line {
            write!(f, " at {}:{line}", self.script)?;
            if let Some(column) = self.column {
                write!(f, ":{column}")?;
            }
        }
        write!(f, ": {}", self.message)
    }
}

impl ScriptRunError {
    fn from_eval(error: &EvalAltResult, script: &str) -> ScriptRunError {
        let position = error.position();
        let (line, column) = (position.line(), position.position());
        let message = error.to_string();
        let (kind, payload) = match error {
            EvalAltResult::ErrorRuntime(value, _) => {
                let payload = serde_json::to_value(value).unwrap_or(Value::Null);
                let kind = payload
                    .get("kind")
                    .and_then(Value::as_str)
                    .unwrap_or("RuntimeError")
                    .to_string();
                (kind, payload)
            }
            EvalAltResult::ErrorParsing(..) => ("SyntaxError".to_string(), Value::Null),
            EvalAltResult::ErrorTerminated(token, _) => (
                "BudgetExceeded".to_string(),
                serde_json::to_value(token).unwrap_or(Value::Null),
            ),
            EvalAltResult::ErrorTooManyOperations(_) => ("BudgetExceeded".to_string(), Value::Null),
            _ => ("ArgumentError".to_string(), Value::Null),
        };
        let message = payload
            .get("message")
            .and_then(Value::as_str)
            .map(str::to_string)
            .unwrap_or(message);
        ScriptRunError {
            kind,
            message,
            script: script.to_string(),
            line,
            column,
            payload: Box::new(payload),
        }
    }
}

/// An engine configured for the popkern dialect.
pub struct ScriptHost {
    engine: Engine,
    printed: Rc<RefCell<Vec<String>>>,
    options: HostOptions,
}

impl ScriptHost {
    /// A host with no budgets.
    pub fn new() -> ScriptHost {
        ScriptHost::with_options(HostOptions::default())
    }

    /// A host with the given budgets.
    pub fn with_options(options: HostOptions) -> ScriptHost {
        let printed: Rc<RefCell<Vec<String>>> = Rc::new(RefCell::new(Vec::new()));
        let sink = Rc::clone(&printed);
        let mut engine = build_engine(move |line| sink.borrow_mut().push(line.to_string()));
        if let Some(max) = options.max_operations {
            engine.set_max_operations(max);
        }
        ScriptHost {
            engine,
            printed,
            options,
        }
    }

    /// The underlying engine (for hosts that add their own functions).
    pub fn engine(&self) -> &Engine {
        &self.engine
    }

    /// Mutable engine access.
    pub fn engine_mut(&mut self) -> &mut Engine {
        &mut self.engine
    }

    /// Parse a script without running it.
    pub fn compile(&self, source: &str, script: &str) -> Result<AST, ScriptRunError> {
        self.engine.compile(source).map_err(|e| {
            let error: EvalAltResult = e.into();
            ScriptRunError::from_eval(&error, script)
        })
    }

    /// Run `source`, named `script` in diagnostics.
    pub fn run(&mut self, source: &str, script: &str) -> Result<RunOutcome, ScriptRunError> {
        let ast = self.compile(source, script)?;
        self.run_ast(&ast, script)
    }

    /// Run a compiled script.
    pub fn run_ast(&mut self, ast: &AST, script: &str) -> Result<RunOutcome, ScriptRunError> {
        self.printed.borrow_mut().clear();
        let watch = Stopwatch::start();
        if let Some(timeout) = self.options.timeout {
            let budget = watch;
            self.engine.on_progress(move |_| {
                if budget.elapsed().is_some_and(|elapsed| elapsed > timeout) {
                    Some(Dynamic::from(format!(
                        "wall-clock budget of {timeout:?} exceeded"
                    )))
                } else {
                    None
                }
            });
        }
        let mut scope = Scope::new();
        let returned = self
            .engine
            .eval_ast_with_scope::<Dynamic>(&mut scope, ast)
            .map_err(|e| ScriptRunError::from_eval(&e, script))?;
        let in_scope: Vec<Dynamic> = scope.iter().map(|(_, _, value)| value).collect();
        drop(scope);
        let model = returned
            .try_cast::<ScriptModel>()
            .or_else(|| {
                in_scope
                    .into_iter()
                    .rev()
                    .find_map(|value| value.try_cast::<ScriptModel>())
            })
            .map(ScriptModel::into_inner);
        Ok(RunOutcome {
            model,
            printed: self.printed.borrow().clone(),
            elapsed: watch.elapsed(),
        })
    }
}

impl Default for ScriptHost {
    fn default() -> Self {
        ScriptHost::new()
    }
}

/// Position helper for tests and hosts.
pub fn position_of(error: &EvalAltResult) -> Position {
    error.position()
}
