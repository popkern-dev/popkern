//! Assembly parameters (PK-638): one value that drives a parameter in
//! several parts — a gearbox's face width in the sun, the planet and
//! the ring.
//!
//! An [`AssemblyParam`] holds an SI value and the unit it is shown in.
//! A [`Drive`] binds it to a declared parameter of one part definition.
//! Setting the assembly parameter sets every driven part parameter (in
//! that parameter's own unit) and replays each driven part once; every
//! instance of those parts follows. Driven part parameters are not
//! edited on their own (a host shows them read-only): the assembly
//! parameter is their single source.

use popkern_history::unit::Unit;
use serde::{Deserialize, Serialize};

use crate::assembly::{Assembly, DefinitionId};
use crate::error::AssemblyError;
use crate::names::validate_name;

/// One assembly-level parameter.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AssemblyParam {
    /// Unique among the assembly's parameters.
    pub name: String,
    /// The value in SI units.
    pub value: f64,
    /// The unit it is shown and set in.
    pub unit: Unit,
}

impl AssemblyParam {
    /// The value in its display unit.
    pub fn display_value(&self) -> f64 {
        self.unit.from_si(self.value)
    }
}

/// An assembly parameter driving a part definition's parameter.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Drive {
    /// The assembly parameter.
    pub param: String,
    /// The driven part definition's name.
    pub part: String,
    /// The part's parameter (its qualified name in that document).
    pub part_param: String,
}

impl Assembly {
    /// The assembly parameters, in declaration order.
    pub fn params(&self) -> &[AssemblyParam] {
        &self.params
    }

    /// The drives, in declaration order.
    pub fn drives(&self) -> &[Drive] {
        &self.drives
    }

    /// Declare an assembly parameter (`value` in `unit`).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::InvalidName`]; [`AssemblyError::Param`] when the
    /// name is taken or the value is not finite.
    pub fn add_param(&mut self, name: &str, value: f64, unit: Unit) -> Result<(), AssemblyError> {
        validate_name(name, "parameter")?;
        if self.params.iter().any(|p| p.name == name) {
            return Err(param_error(
                name,
                "an assembly parameter with this name already exists",
            ));
        }
        if !value.is_finite() {
            return Err(param_error(name, "the value must be finite"));
        }
        self.params.push(AssemblyParam {
            name: name.to_string(),
            value: unit.to_si(value),
            unit,
        });
        Ok(())
    }

    /// Bind assembly parameter `param` to `part_param` of the part
    /// `definition`, and apply its value to the part now.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Param`] for an unknown assembly parameter, a part
    /// with no such parameter, or a part parameter already driven;
    /// [`AssemblyError::WrongDefinitionKind`] for a sub-assembly;
    /// [`AssemblyError::Script`] when the part does not replay.
    pub fn drive(
        &mut self,
        param: &str,
        definition: DefinitionId,
        part_param: &str,
    ) -> Result<(), AssemblyError> {
        let Some(value) = self
            .params
            .iter()
            .find(|p| p.name == param)
            .map(|p| p.value)
        else {
            return Err(param_error(param, "no assembly parameter with this name"));
        };
        let part = self.definition(definition)?.name.clone();
        if self
            .drives
            .iter()
            .any(|d| d.part == part && d.part_param == part_param)
        {
            return Err(param_error(
                param,
                &format!("'{part}:{part_param}' is already driven"),
            ));
        }
        apply(self, definition, &part, part_param, value)?;
        self.drives.push(Drive {
            param: param.to_string(),
            part,
            part_param: part_param.to_string(),
        });
        Ok(())
    }

    /// Set assembly parameter `name` to `value` (in its unit): every
    /// driven part parameter follows and each driven part replays once.
    /// Returns the parts rebuilt.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Param`] for an unknown parameter or a non-finite
    /// value; [`AssemblyError::Script`] when a driven part fails.
    pub fn set_assembly_param(
        &mut self,
        name: &str,
        value: f64,
    ) -> Result<Vec<String>, AssemblyError> {
        if !value.is_finite() {
            return Err(param_error(name, "the value must be finite"));
        }
        let Some(index) = self.params.iter().position(|p| p.name == name) else {
            return Err(param_error(name, "no assembly parameter with this name"));
        };
        let si = self.params[index].unit.to_si(value);
        self.params[index].value = si;
        let drives: Vec<Drive> = self
            .drives
            .iter()
            .filter(|d| d.param == name)
            .cloned()
            .collect();
        let mut rebuilt = Vec::new();
        for drive in drives {
            let Some(id) = self.definition_named(&drive.part) else {
                return Err(param_error(
                    name,
                    &format!("driven part '{}' is gone", drive.part),
                ));
            };
            apply(self, id, &drive.part, &drive.part_param, si)?;
            if !rebuilt.contains(&drive.part) {
                rebuilt.push(drive.part);
            }
        }
        Ok(rebuilt)
    }

    /// Whether `part:part_param` is driven by an assembly parameter
    /// (and by which).
    pub fn driven_by(&self, part: &str, part_param: &str) -> Option<&str> {
        self.drives
            .iter()
            .find(|d| d.part == part && d.part_param == part_param)
            .map(|d| d.param.as_str())
    }
}

/// Set the part's parameter to the SI value (converted to the part
/// parameter's own unit) and replay the part.
fn apply(
    assembly: &mut Assembly,
    id: DefinitionId,
    part: &str,
    part_param: &str,
    si: f64,
) -> Result<(), AssemblyError> {
    let script = |source| AssemblyError::Script {
        definition: part.to_string(),
        source,
    };
    let model = assembly.part_model_mut(id)?;
    let Some(unit) = model
        .params()
        .into_iter()
        .find(|p| p.name == part_param)
        .map(|p| p.unit)
    else {
        return Err(param_error(
            part_param,
            &format!("part '{part}' declares no parameter '{part_param}'"),
        ));
    };
    model
        .set_param(part_param, unit.from_si(si))
        .map_err(script)?;
    model.replay().map_err(script)?;
    Ok(())
}

fn param_error(name: &str, detail: &str) -> AssemblyError {
    AssemblyError::Param {
        name: name.to_string(),
        detail: detail.to_string(),
    }
}
