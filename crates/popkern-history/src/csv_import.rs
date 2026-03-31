//! CSV import implementation for configurations.
//!
//! This module is not public API — it is a submodule of [`crate::csv`].
//! All pub items here are re-exported through `crate::csv`.

use std::collections::HashMap;

use crate::FeatureTree;
use crate::config::{ConfigId, ParamOverride};
use crate::csv::{CsvImportResult, parse_csv_row, parse_literal};
use crate::error::HistoryError;
use crate::node::FeatureId;
use crate::parameter::ParameterId;
use crate::unit::Unit;

// ---------------------------------------------------------------------------
// Column-kind discriminant used during import
// ---------------------------------------------------------------------------

/// Column kind resolved during CSV header parsing.
#[derive(Debug)]
pub(super) enum ColKind {
    /// `name` metadata column.
    Name,
    /// `description` metadata column.
    Description,
    /// `parent` metadata column.
    Parent,
    /// A parameter override column (resolved to its [`ParameterId`]).
    Param(ParameterId),
    /// A feature suppression column (resolved to its [`FeatureId`]).
    Suppression(FeatureId),
}

// ---------------------------------------------------------------------------
// Import entry point
// ---------------------------------------------------------------------------

impl FeatureTree {
    /// Import configurations from a CSV string.
    ///
    /// Rows with the name `Default` update base parameter values rather than
    /// creating a configuration override.  All other rows create or update named
    /// configurations (matched by name).
    ///
    /// Non-fatal issues (unknown columns, bad cells) are collected in
    /// [`CsvImportResult::warnings`] rather than aborting the import.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError`] only for structural failures such as
    /// duplicate config names; individual cell errors are recorded as warnings.
    pub fn import_configs_csv(&mut self, csv: &str) -> Result<CsvImportResult, HistoryError> {
        let mut result = CsvImportResult {
            created: Vec::new(),
            updated: Vec::new(),
            warnings: Vec::new(),
        };

        // Capture undo cursor so we can wrap all import operations into a
        // single compound undo step.
        let undo_start = self.undo_cursor;

        // Split into lines, handling both \r\n and \n.
        let lines: Vec<&str> = csv.lines().map(|l| l.trim_end_matches('\r')).collect();

        if lines.is_empty() {
            return Ok(result);
        }

        // --- Parse header ---
        let header_fields = parse_csv_row(lines[0]);
        let col_kinds = self.resolve_header(&header_fields, &mut result.warnings);

        // Build a name -> ConfigId map for the currently existing configs.
        let mut name_to_id: HashMap<String, ConfigId> = self
            .configs()
            .into_iter()
            .filter_map(|id| self.config(id).map(|c| (c.name.clone(), id)))
            .collect();

        // Deferred parent assignments: all configs must be created before
        // parents can be resolved.
        struct PendingParent {
            config_id: ConfigId,
            parent_name: String,
        }
        let mut pending_parents: Vec<PendingParent> = Vec::new();

        // First pass: create / identify configs and apply overrides.
        for (row_idx, line) in lines[1..].iter().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let fields = parse_csv_row(line);

            let name = col_kinds
                .iter()
                .zip(fields.iter())
                .find(|(k, _)| matches!(k, ColKind::Name))
                .map(|(_, v)| v.trim().to_string())
                .unwrap_or_default();

            if name.is_empty() {
                result
                    .warnings
                    .push(format!("row {}: empty name, skipping", row_idx + 2));
                continue;
            }

            let description = col_kinds
                .iter()
                .zip(fields.iter())
                .find(|(k, _)| matches!(k, ColKind::Description))
                .map(|(_, v)| v.trim().to_string())
                .unwrap_or_default();

            let parent_name = col_kinds
                .iter()
                .zip(fields.iter())
                .find(|(k, _)| matches!(k, ColKind::Parent))
                .map(|(_, v)| v.trim().to_string())
                .unwrap_or_default();

            if name == "Default" {
                self.apply_default_row(&col_kinds, &fields, row_idx + 2, &mut result.warnings)?;
                continue;
            }

            let cfg_id = if let Some(&existing_id) = name_to_id.get(&name) {
                result.updated.push(existing_id);
                existing_id
            } else {
                let new_id = self.add_config(&name, &description)?;
                name_to_id.insert(name.clone(), new_id);
                result.created.push(new_id);
                new_id
            };

            if !parent_name.is_empty() {
                pending_parents.push(PendingParent {
                    config_id: cfg_id,
                    parent_name,
                });
            }

            self.apply_config_row(
                cfg_id,
                &col_kinds,
                &fields,
                row_idx + 2,
                &mut result.warnings,
            )?;
        }

        // Second pass: assign parents (all configs now exist).
        for pp in pending_parents {
            if let Some(&parent_id) = name_to_id.get(&pp.parent_name) {
                if let Some(cfg) = self.configs.get_mut(pp.config_id) {
                    cfg.parent = Some(parent_id);
                }
            } else {
                result.warnings.push(format!(
                    "parent '{}' not found for config '{}'",
                    pp.parent_name,
                    self.config(pp.config_id)
                        .map(|c| c.name.as_str())
                        .unwrap_or("?")
                ));
            }
        }

        // Wrap all individual undo entries from this import into a single
        // compound command so the entire import is one undo step.
        if self.undo_cursor > undo_start {
            let individual_cmds: Vec<_> =
                self.undo_log.drain(undo_start..self.undo_cursor).collect();
            self.undo_log
                .push(crate::undo::TreeCommand::Compound(individual_cmds));
            self.undo_cursor = self.undo_log.len();
        }

        Ok(result)
    }

    /// Resolve header fields into [`ColKind`] variants.
    ///
    /// Unknown columns produce a warning and fall back to [`ColKind::Name`]
    /// (harmless sentinel — the processing loops skip duplicate Name entries).
    pub(super) fn resolve_header(
        &self,
        fields: &[String],
        warnings: &mut Vec<String>,
    ) -> Vec<ColKind> {
        let label_to_fid: HashMap<String, FeatureId> = self
            .timeline()
            .iter()
            .filter_map(|&fid| self.feature(fid).map(|n| (n.label.clone(), fid)))
            .collect();

        fields
            .iter()
            .map(|h| {
                let h = h.trim();
                match h {
                    "name" => ColKind::Name,
                    "description" => ColKind::Description,
                    "parent" => ColKind::Parent,
                    _ => {
                        if let Some(label) = h.strip_suffix(".suppressed") {
                            if let Some(&fid) = label_to_fid.get(label) {
                                return ColKind::Suppression(fid);
                            }
                            warnings.push(format!(
                                "unknown suppression column '{h}' (feature '{label}' not found)"
                            ));
                        }
                        if let Some(pid) = self.parameter_table.resolve_name(h) {
                            return ColKind::Param(pid);
                        }
                        warnings.push(format!("unknown column '{h}', will be skipped"));
                        ColKind::Name
                    }
                }
            })
            .collect()
    }

    /// Apply a `Default` row: update base parameter values.
    fn apply_default_row(
        &mut self,
        col_kinds: &[ColKind],
        fields: &[String],
        row_num: usize,
        warnings: &mut Vec<String>,
    ) -> Result<(), HistoryError> {
        for (kind, cell) in col_kinds.iter().zip(fields.iter()) {
            let cell = cell.trim();
            if cell.is_empty() {
                continue;
            }
            if let ColKind::Param(pid) = kind {
                if cell.starts_with('=') {
                    warnings.push(format!(
                        "row {row_num}: Default row param cell '{cell}' is an expression; \
                         only literals allowed for base values, skipping"
                    ));
                    continue;
                }
                match parse_literal(cell) {
                    Some(si_val) => {
                        let unit = self
                            .parameter_table
                            .get(*pid)
                            .map(|p| p.unit)
                            .unwrap_or(Unit::Ratio);
                        if let Some(p) = self.parameter_table.get_mut(*pid)
                            && !p.read_only
                        {
                            p.value = si_val;
                            p.unit = unit;
                        }
                    }
                    None => {
                        warnings.push(format!(
                            "row {row_num}: cannot parse literal '{cell}', skipping"
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// Apply parameter and suppression overrides for a named config row.
    fn apply_config_row(
        &mut self,
        cfg_id: ConfigId,
        col_kinds: &[ColKind],
        fields: &[String],
        row_num: usize,
        warnings: &mut Vec<String>,
    ) -> Result<(), HistoryError> {
        for (kind, cell) in col_kinds.iter().zip(fields.iter()) {
            let cell = cell.trim();
            if cell.is_empty() {
                continue;
            }
            match kind {
                ColKind::Param(pid) => {
                    let pid = *pid;
                    let ov = if let Some(expr_str) = cell.strip_prefix('=') {
                        match crate::parse::parse_expression(expr_str, &self.parameter_table) {
                            Ok(expr) => Some(ParamOverride::Expression(expr)),
                            Err(e) => {
                                warnings.push(format!(
                                    "row {row_num}: expression parse error in '{cell}': {e}"
                                ));
                                None
                            }
                        }
                    } else {
                        match parse_literal(cell) {
                            Some(si_val) => Some(ParamOverride::Literal(si_val)),
                            None => {
                                warnings.push(format!(
                                    "row {row_num}: cannot parse literal '{cell}', skipping"
                                ));
                                None
                            }
                        }
                    };
                    if let Some(override_val) = ov {
                        self.set_config_override(cfg_id, pid, override_val)?;
                    }
                }
                ColKind::Suppression(fid) => {
                    let fid = *fid;
                    match cell {
                        "true" => self.set_config_suppression(cfg_id, fid, true),
                        "false" => self.set_config_suppression(cfg_id, fid, false),
                        _ => {
                            warnings.push(format!(
                                "row {row_num}: invalid suppression value '{cell}', \
                                 expected 'true' or 'false'"
                            ));
                        }
                    }
                }
                _ => {} // metadata or sentinel columns — already handled above
            }
        }
        Ok(())
    }
}
