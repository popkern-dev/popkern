//! CSV import and export for named configurations.
//!
//! Provides round-trip serialisation of the configuration design table so that
//! users can edit parameter overrides in a spreadsheet application and import
//! them back.
//!
//! # Format
//!
//! ```text
//! name,description,parent,width,fillet1.radius,Chamfer1.suppressed
//! Default,Base design,,50mm,=width * 0.1,false
//! Small,Residential,,25mm,=width * 0.1,true
//! Large,Commercial,,100mm,=width * 0.05,false
//! ```
//!
//! The first three columns are always `name`, `description`, and `parent`.
//! Subsequent columns are parameter overrides (sorted by qualified name) and
//! then suppression overrides (one column per feature that has any suppression
//! override in any configuration, named `<feature_label>.suppressed`).
//!
//! RFC 4180 quoting is used: cells containing `,`, `"`, or newlines are
//! wrapped in double quotes; internal double quotes become `""`.

#[path = "csv_import.rs"]
mod csv_import;

use std::collections::BTreeMap;

use crate::FeatureTree;
use crate::config::ParamOverride;
use crate::node::FeatureId;
use crate::parameter::ParameterId;
use crate::unit::Unit;

// ---------------------------------------------------------------------------
// Public result type
// ---------------------------------------------------------------------------

/// Result of importing configurations from CSV.
#[derive(Debug)]
pub struct CsvImportResult {
    /// Configurations that were created during the import.
    pub created: Vec<crate::config::ConfigId>,
    /// Configurations that were updated (matched by name).
    pub updated: Vec<crate::config::ConfigId>,
    /// Non-fatal warnings (unparseable cells, unknown columns, etc.).
    pub warnings: Vec<String>,
}

// ---------------------------------------------------------------------------
// CSV primitive helpers (pub(crate) so csv_import can use them)
// ---------------------------------------------------------------------------

/// Wrap a cell value in RFC 4180 quotes when necessary.
pub(crate) fn csv_quote(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Parse a single CSV row, handling RFC 4180 quoting and `""` escape.
pub(crate) fn parse_csv_row(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();

    while let Some(ch) = chars.next() {
        if in_quotes {
            if ch == '"' {
                if chars.peek() == Some(&'"') {
                    chars.next(); // consume second quote
                    field.push('"');
                } else {
                    in_quotes = false;
                }
            } else {
                field.push(ch);
            }
        } else if ch == '"' {
            in_quotes = true;
        } else if ch == ',' {
            fields.push(field.clone());
            field.clear();
        } else {
            field.push(ch);
        }
    }
    fields.push(field);
    fields
}

/// Parse a display-value string such as `"50mm"`, `"45deg"`, or `"3.5"` into a
/// value in SI base units.  Returns `None` on failure.
pub(crate) fn parse_literal(s: &str) -> Option<f64> {
    let s = s.trim();
    let num_end = s
        .find(|c: char| c.is_alphabetic() || c == '_')
        .unwrap_or(s.len());
    let num_str = s[..num_end].trim();
    let unit_str = s[num_end..].trim();

    let raw: f64 = num_str.parse().ok()?;

    // Match unit suffixes (longest first to avoid "m" matching "mm").
    let unit = match unit_str {
        "mm" => Unit::Millimeters,
        "cm" => Unit::Centimeters,
        "deg" => Unit::Degrees,
        "rad" => Unit::Radians,
        "in" => Unit::Inches,
        "m" => Unit::Meters,
        "" => Unit::Ratio,
        _ => return None,
    };

    Some(unit.to_si(raw))
}

/// Format a parameter value for export: convert from SI, then append unit suffix.
fn format_param_value(si_value: f64, unit: Unit) -> String {
    let display = unit.from_si(si_value);
    let suffix = unit.suffix();
    let s = format!("{display:.10}");
    let s = if s.contains('.') {
        let s = s.trim_end_matches('0');
        s.trim_end_matches('.')
    } else {
        s.as_str()
    };
    if suffix.is_empty() {
        s.to_string()
    } else {
        format!("{s}{suffix}")
    }
}

// ---------------------------------------------------------------------------
// Export
// ---------------------------------------------------------------------------

impl FeatureTree {
    /// Export all configurations as a CSV string.
    ///
    /// The CSV uses RFC 4180 quoting.  The first row is a header; the second
    /// row is always `Default` (showing base parameter values); subsequent rows
    /// are named configurations.
    ///
    /// Column order: `name`, `description`, `parent`, then parameter columns
    /// sorted alphabetically by qualified name, then suppression columns sorted
    /// alphabetically by feature label.
    pub fn export_configs_csv(&self) -> String {
        let table = &self.parameter_table;

        // Collect all non-read-only parameters, sorted by qualified name.
        let mut param_cols: Vec<(String, ParameterId)> = table
            .all_params()
            .filter(|(_, p)| !p.read_only)
            .filter_map(|(id, _)| table.qualified_name(id).map(|name| (name, id)))
            .collect();
        param_cols.sort_by(|a, b| a.0.cmp(&b.0));

        // Collect features that appear in any suppression map.
        let mut suppressed_features: BTreeMap<String, FeatureId> = BTreeMap::new();
        for cfg_id in self.configs() {
            if let Some(suppressions) = self.config_suppressions(cfg_id) {
                for &fid in suppressions.keys() {
                    if let Some(node) = self.feature(fid) {
                        suppressed_features.entry(node.label.clone()).or_insert(fid);
                    }
                }
            }
        }

        // Build header row.
        let mut headers: Vec<String> = vec![
            "name".to_string(),
            "description".to_string(),
            "parent".to_string(),
        ];
        for (name, _) in &param_cols {
            headers.push(name.clone());
        }
        for label in suppressed_features.keys() {
            headers.push(format!("{label}.suppressed"));
        }
        let header_row = headers
            .iter()
            .map(|h| csv_quote(h))
            .collect::<Vec<_>>()
            .join(",");

        // Build the Default row with base parameter values.
        let mut rows = vec![header_row];
        let mut default_cells: Vec<String> =
            vec!["Default".to_string(), String::new(), String::new()];
        for (_, pid) in &param_cols {
            if let Some(param) = table.get(*pid) {
                default_cells.push(format_param_value(param.value, param.unit));
            } else {
                default_cells.push(String::new());
            }
        }
        for _ in suppressed_features.keys() {
            default_cells.push(String::new());
        }
        rows.push(
            default_cells
                .iter()
                .map(|c| csv_quote(c))
                .collect::<Vec<_>>()
                .join(","),
        );

        // Build one row per config.
        for cfg_id in self.configs() {
            let cfg = match self.config(cfg_id) {
                Some(c) => c,
                None => continue,
            };

            let parent_name = cfg
                .parent
                .and_then(|pid| self.config(pid))
                .map(|p| p.name.clone())
                .unwrap_or_default();

            let overrides = cfg.overrides.clone();
            let suppressions = cfg.suppression.clone();

            let mut cells: Vec<String> =
                vec![cfg.name.clone(), cfg.description.clone(), parent_name];

            for (_, pid) in &param_cols {
                let cell = if let Some(ov) = overrides.get(pid) {
                    match ov {
                        ParamOverride::Literal(v) => {
                            let unit = table.get(*pid).map(|p| p.unit).unwrap_or(Unit::Ratio);
                            format_param_value(*v, unit)
                        }
                        ParamOverride::Expression(expr) => {
                            format!("={}", crate::format::format_expression(expr, table))
                        }
                    }
                } else {
                    String::new()
                };
                cells.push(cell);
            }

            for fid in suppressed_features.values() {
                let cell = if let Some(&suppressed) = suppressions.get(fid) {
                    if suppressed {
                        "true".to_string()
                    } else {
                        "false".to_string()
                    }
                } else {
                    String::new()
                };
                cells.push(cell);
            }

            rows.push(
                cells
                    .iter()
                    .map(|c| csv_quote(c))
                    .collect::<Vec<_>>()
                    .join(","),
            );
        }

        rows.join("\n")
    }
}

#[cfg(test)]
#[path = "csv_tests.rs"]
mod tests;
