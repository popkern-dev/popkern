//! `fillet_each` on a script body (PK-602), split from `body_api` for
//! the file budget: parse `[[edge selector, radius], …]` groups and
//! record one grouped fillet.

use popkern_script::EdgeSelector;
use rhai::{Dynamic, EvalAltResult};

use super::ScriptBody;
use crate::args;
use crate::error::IntoScriptResult;
use crate::recipe_api::fillet_options;

impl ScriptBody {
    pub(super) fn fillet_each(
        &self,
        groups: rhai::Array,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let mut parsed = Vec::with_capacity(groups.len());
        for (i, group) in groups.into_iter().enumerate() {
            let shape_error = || -> Box<EvalAltResult> {
                format!("fillet_each: group {i} must be [edge selector, radius]").into()
            };
            let pair = group.try_cast::<rhai::Array>().ok_or_else(shape_error)?;
            if pair.len() != 2 {
                return Err(shape_error());
            }
            let selector = pair[0]
                .clone()
                .try_cast::<EdgeSelector>()
                .ok_or_else(shape_error)?;
            let radius = args::scalar("radius", &pair[1])?;
            parsed.push((selector, radius));
        }
        let (recipe, law, label) = fillet_options(options)?;
        if law != popkern_script::RadiusLaw::Constant {
            return Err(
                "fillet_each: radius laws are per fillet; fillet_each groups take constant radii"
                    .into(),
            );
        }
        let body = self
            .body
            .fillet_each(&mut self.model.borrow_mut(), parsed, recipe, &label)
            .script()?;
        Ok(self.served(body))
    }
}
