//! Recorder errors as script values. Every `ScriptError` crosses into
//! Rhai as a runtime error whose value is the spec-§10 payload map
//! (`ScriptError::payload`, shared by every skin) with `kind`,
//! `message`, and the class's fields; the engine stamps the script
//! position of the call that raised it, and `try`/`catch` receives the
//! same map.

use popkern_script::ScriptError;
use rhai::{Dynamic, EvalAltResult, Map, Position};
use serde_json::{Value, json};

/// Build a script runtime error from a §10 payload.
pub fn raise(payload: Value) -> Box<EvalAltResult> {
    let value =
        rhai::serde::to_dynamic(&payload).unwrap_or_else(|_| Dynamic::from(payload.to_string()));
    Box::new(EvalAltResult::ErrorRuntime(value, Position::NONE))
}

/// A skin-level `ArgumentError` (raised before the kernel runs).
pub fn argument_error(argument: &str, reason: impl Into<String>) -> Box<EvalAltResult> {
    let reason = reason.into();
    raise(json!({
        "kind": "ArgumentError",
        "message": format!("argument '{argument}': {reason}"),
        "argument": argument,
        "reason": reason,
    }))
}

/// The §10 payload of a recorder error (the shared projection).
pub fn payload(error: &ScriptError) -> Value {
    error.payload()
}

/// Convert a recorder error into the script runtime error.
pub fn script_error(error: ScriptError) -> Box<EvalAltResult> {
    raise(error.payload())
}

/// `Result` adaptor: `recorder_call().script()?`.
pub trait IntoScriptResult<T> {
    /// Map the recorder error into a script runtime error.
    fn script(self) -> Result<T, Box<EvalAltResult>>;
}

impl<T> IntoScriptResult<T> for Result<T, ScriptError> {
    fn script(self) -> Result<T, Box<EvalAltResult>> {
        self.map_err(script_error)
    }
}

/// The raising call's error, unwrapped from the function-call frames a
/// closure adds (`m.tool(|m| …)`, `a.part(…, |m| …)`): the class,
/// payload and position are the call that failed, not the closure call.
pub fn innermost(error: Box<EvalAltResult>) -> Box<EvalAltResult> {
    let mut error = error;
    while let EvalAltResult::ErrorInFunctionCall(_, _, inner, _) = *error {
        error = inner;
    }
    error
}

/// Read the error map back out of a runtime error (host side).
pub fn error_map(error: &EvalAltResult) -> Option<Map> {
    match error {
        EvalAltResult::ErrorRuntime(value, _) => value.clone().try_cast::<Map>(),
        _ => None,
    }
}

/// `Lateral/1`-style text for a serialised `EntityName`.
fn entity_label(name: Option<&Value>) -> Option<String> {
    let role = name?.get("role")?;
    let category = role.get("category")?.as_str()?;
    let index = role.get("index")?.as_u64()?;
    Some(format!("{category}/{index}"))
}

/// Remediation lines a payload carries: the whole-error hint and the
/// per-edge hints of a blend decline, as `action parameter = value:
/// rationale`.
pub fn remediation_lines(payload: &Value) -> Vec<String> {
    fn line(prefix: &str, remediation: &Value) -> Option<String> {
        let action = remediation.get("action")?.as_str()?;
        let parameter = remediation
            .get("parameter")
            .and_then(Value::as_str)
            .unwrap_or("");
        let value = remediation.get("value").and_then(Value::as_f64);
        let rationale = remediation
            .get("rationale")
            .and_then(Value::as_str)
            .unwrap_or("");
        let mut text = format!("{prefix}{action}");
        if !parameter.is_empty() {
            text.push_str(&format!(" {parameter}"));
        }
        if let Some(v) = value {
            text.push_str(&format!(" = {v}"));
        }
        if !rationale.is_empty() {
            text.push_str(&format!(": {rationale}"));
        }
        Some(text)
    }
    let error = payload.get("error").unwrap_or(payload);
    let mut lines = Vec::new();
    if let Some(r) = error.get("remediation").and_then(|r| line("", r)) {
        lines.push(r);
    }
    if let Some(failures) = error.get("failures").and_then(Value::as_array) {
        for failure in failures {
            let edge = entity_label(failure.get("edge_name")).unwrap_or_else(|| {
                failure
                    .get("edge")
                    .and_then(Value::as_str)
                    .unwrap_or("edge")
                    .to_string()
            });
            if let Some(r) = failure
                .get("remediation")
                .and_then(|r| line(&format!("{edge}: "), r))
            {
                lines.push(r);
            }
        }
    }
    lines
}
