//! Parameter expressions and configurations (spec §4.2–§4.3): thin,
//! name-keyed wrappers over `popkern_history::config_api`.

use popkern_history::config::{ConfigId, ParamOverride};
use popkern_history::unit::Unit;

use crate::error::ScriptError;
use crate::model::{FeatureHandle, Model};

impl Model {
    /// Declare a global parameter driven by an expression.
    pub fn param_expr(&mut self, name: &str, text: &str, unit: Unit) -> Result<(), ScriptError> {
        let expression = popkern_history::parse::parse_expression(text, self.tree.parameters())
            .map_err(|e| ScriptError::Argument {
                argument: name.to_string(),
                reason: format!("expression '{text}': {e}"),
            })?;
        let id = self.tree.add_global_parameter(name, 0.0, unit)?;
        self.tree.set_parameter_expression(id, expression)?;
        Ok(())
    }

    fn config_id(&self, name: &str) -> Result<ConfigId, ScriptError> {
        self.tree
            .configs()
            .into_iter()
            .find(|&id| self.tree.config(id).map(|c| c.name.as_str()) == Some(name))
            .ok_or_else(|| popkern_history::HistoryError::ConfigNotFound.into())
    }

    fn param_id(&self, name: &str) -> Result<popkern_history::parameter::ParameterId, ScriptError> {
        self.tree
            .parameters()
            .resolve_name(name)
            .ok_or_else(|| ScriptError::Argument {
                argument: "parameter".into(),
                reason: format!("unknown parameter '{name}'"),
            })
    }

    /// Add a configuration.
    pub fn add_config(&mut self, name: &str, description: &str) -> Result<(), ScriptError> {
        self.tree.add_config(name, description)?;
        Ok(())
    }

    /// Activate a configuration (`None` = base) and replay.
    pub fn set_active_config(&mut self, name: Option<&str>) -> Result<(), ScriptError> {
        let id = match name {
            Some(n) => Some(self.config_id(n)?),
            None => None,
        };
        self.tree.set_active_config(id);
        self.tree.replay(self.tolerance)?;
        Ok(())
    }

    /// The active configuration's name.
    pub fn active_config(&self) -> Option<String> {
        self.tree
            .active_config()
            .and_then(|id| self.tree.config(id))
            .map(|c| c.name.clone())
    }

    /// Override a parameter (display value) in a configuration.
    pub fn set_override(
        &mut self,
        config: &str,
        param: &str,
        value: f64,
    ) -> Result<(), ScriptError> {
        let cid = self.config_id(config)?;
        let pid = self.param_id(param)?;
        let unit = self
            .tree
            .parameters()
            .get(pid)
            .map(|p| p.unit)
            .unwrap_or(self.unit);
        self.tree
            .set_config_override(cid, pid, ParamOverride::Literal(unit.to_si(value)))?;
        Ok(())
    }

    /// Remove a configuration override.
    pub fn clear_override(&mut self, config: &str, param: &str) -> Result<(), ScriptError> {
        let cid = self.config_id(config)?;
        let pid = self.param_id(param)?;
        self.tree.clear_config_override(cid, pid);
        Ok(())
    }

    /// Suppress (or unsuppress) a feature in a configuration.
    pub fn set_suppression(
        &mut self,
        config: &str,
        feature: &FeatureHandle,
        suppressed: bool,
    ) -> Result<(), ScriptError> {
        let cid = self.config_id(config)?;
        self.tree
            .set_config_suppression(cid, feature.id, suppressed);
        Ok(())
    }

    /// Configuration names.
    pub fn configs(&self) -> Vec<String> {
        self.tree
            .configs()
            .into_iter()
            .filter_map(|id| self.tree.config(id).map(|c| c.name.clone()))
            .collect()
    }
}
